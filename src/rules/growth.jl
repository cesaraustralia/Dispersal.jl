# Mixins

@premix @columns struct InstrinsicGrowthRate{GR}
    # Field           | Def | Flatn | Limits      | Description
    intrinsicrate::GR | 0.1 | true  | (0.0, 10.0) | "Intrinsic rate of growth"
end

@premix @columns struct CarryCap{CC}
    carrycap::CC | 100000   | true  | (0.0, 1e9) | "Carrying capacity for each cell. Not currently scaled by area."
end

@premix @columns struct Layers{L}
    layers::L    | ()       | false | _ | "Additional data layers"
end


# Type declarations

"""
Extends AbstractCellRule for rules of growth dynamics

For best performance these should be chained with other 
AbstractCellRule or following an AbstractNeighborhoodRule.
"""
abstract type AbstractGrowthRule <: AbstractCellRule end


# Euler method solvers

"""
Simple fixed exponential growth rate solved with Euler method
$(FIELDDOCTABLE)
"""
@InstrinsicGrowthRate struct EulerExponentialGrowth{} <: AbstractGrowthRule end

"""
Simple fixed logistic growth rate solved with Euler method
$(FIELDDOCTABLE)
"""
@CarryCap @InstrinsicGrowthRate struct EulerLogisticGrowth{} <: AbstractGrowthRule end

"""
Exponential growth based on a suitability layer solved with Euler method
$(FIELDDOCTABLE)
"""
@Layers struct SuitabilityEulerExponentialGrowth{} <: AbstractGrowthRule end

"""
Logistic growth based on a suitability layer solved with Euler method
$(FIELDDOCTABLE)
"""
@CarryCap @Layers struct SuitabilityEulerLogisticGrowth{} <: AbstractGrowthRule end


# Exact growth solutions

"""
Simple fixed exponential growth rate using exact solution
$(FIELDDOCTABLE)
"""
@InstrinsicGrowthRate struct ExactExponentialGrowth{} <: AbstractGrowthRule end

"""
Simple fixed logistic growth rate using exact solution
$(FIELDDOCTABLE)
"""
@CarryCap @InstrinsicGrowthRate struct ExactLogisticGrowth{} <: AbstractGrowthRule end

"""
Exponential growth based on a suitability layer using exact solution
$(FIELDDOCTABLE)
"""
@Layers struct SuitabilityExactExponentialGrowth{} <: AbstractGrowthRule end

"""
Logistic growth based on a suitability layer using exact solution
$(FIELDDOCTABLE)
"""
@CarryCap @Layers struct SuitabilityExactLogisticGrowth{} <: AbstractGrowthRule end

"""
Simple suitability layer mask
$(FIELDDOCTABLE)
"""
@Layers struct SuitabilityMask{ST} <: AbstractGrowthRule
    threshold::ST | 0.7 | true  | (0.0, 1.0) | "Minimum habitat suitability index."
end


# Rules

# Euler solver rules
@inline applyrule(rule::EulerExponentialGrowth, data, state, args...) = state + state * rule.intrinsicrate * data.timestep

@inline applyrule(rule::SuitabilityEulerExponentialGrowth, data, state, index, args...) = begin
    state == zero(state) && return state
    intrinsicrate = get_layers(rule, data, index)
    state + intrinsicrate * state * data.timestep # dN = rN * dT
    # max(min(state * intrinsicrate, rule.max), rule.min)
end

# TODO: fix and test logistic growth
@inline applyrule(rule::EulerLogisticGrowth, data, state, args...) =
    state + state * rule.intrinsicrate * (oneunit(state) - state / rule.carrycap) * data.timestep # dN = (1-N/K)rN dT

@inline applyrule(rule::SuitabilityEulerLogisticGrowth, data, state, index, args...) = begin
    state == zero(state) && return state
    intrinsicrate = get_layers(rule, data, index)
    saturation = intrinsicrate > zero(intrinsicrate) ? (oneunit(state) - state / rule.carrycap) : oneunit(state)
    state + state * saturation * intrinsicrate * data.timestep
end


# Exact solution rules

@inline applyrule(rule::ExactExponentialGrowth, data, state, args...) = state * exp(rule.intrinsicrate * data.timestep)

@inline applyrule(rule::SuitabilityExactExponentialGrowth, data, state, index, args...) = begin
    state == zero(state) && return state
    intrinsicrate = get_layers(rule, data, index)
    state * exp(intrinsicrate * data.timestep)
    # max(min(state * intrinsicrate, rule.max), rule.min)
end

@inline applyrule(rule::ExactLogisticGrowth, data, state, args...) =
    (state * rule.carrycap) /
    (state + (rule.carrycap - state) * exp(-rule.intrinsicrate * data.timestep))

@inline applyrule(rule::SuitabilityExactLogisticGrowth, data, state, index, args...) = begin
    state == zero(state) && return state
    intrinsicrate = get_layers(rule, data, index)
    # Saturation only applies with positive growth
    if intrinsicrate > zero(intrinsicrate)
        (state * rule.carrycap) /
        (state + (rule.carrycap - state) * exp(-intrinsicrate * data.timestep))
    else
        state * exp(intrinsicrate * data.timestep)
    end
end

@inline applyrule(rule::SuitabilityMask, data, state, index, args...) = begin
    state == zero(state) && return zero(state)
    get_layers(rule, data, index) >= rule.threshold ? state : zero(state)
end

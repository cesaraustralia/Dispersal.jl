# Mixins

@premix @columns struct InstrinsicGrowthRate{GR}
    # Field      | Defaule  | Flatn | Limits      | Description
    intrinsicrate::GR | 0.1 | true  | (0.0, 10.0) | "Intrinsic rate of growth"
end

@premix @columns struct CarryCap{CC}
    carrycap::CC | 100000.0 | true  | (0.0, 1e9)  | "Carrying capacity for each cell. Not currently scaled by area."
end

@premix @columns struct Layers{L,TI}
    layer::L       | nothing          | false | _ | "Data layer"
    timeinterp::TI | Ref((1, 1, 1.0)) | false | _ | "Precaclulated interpolation indices"
end


# Type declarations

"""
Extends AbstractCellRule for rules of growth dynamics

For best performance these should be chained with other
AbstractCellRule or following an AbstractNeighborhoodRule.
"""
abstract type AbstractGrowthRule <: AbstractCellRule end

abstract type AbstractSuitabilityGrowthRule <: AbstractGrowthRule end

layer(rule::AbstractSuitabilityGrowthRule) = rule.layer 
timeinterp(rule::AbstractSuitabilityGrowthRule) = rule.timeinterp[]

DynamicGrids.precalcrule!(rule::AbstractSuitabilityGrowthRule, data) = 
    rule.timeinterp[] = precalc_time_interpolation(layer(rule), rule, data)


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
@Layers struct SuitabilityEulerExponentialGrowth{} <: AbstractSuitabilityGrowthRule end

"""
Logistic growth based on a suitability layer solved with Euler method
$(FIELDDOCTABLE)
"""
@CarryCap @Layers struct SuitabilityEulerLogisticGrowth{} <: AbstractSuitabilityGrowthRule end


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
@Layers struct SuitabilityExactExponentialGrowth{} <: AbstractSuitabilityGrowthRule end

"""
Logistic growth based on a suitability layer using exact solution
$(FIELDDOCTABLE)
"""
@CarryCap @Layers struct SuitabilityExactLogisticGrowth{} <: AbstractSuitabilityGrowthRule end

"""
Simple suitability layer mask
$(FIELDDOCTABLE)
"""
@Layers struct SuitabilityMask{ST} <: AbstractSuitabilityGrowthRule
    threshold::ST | 0.7 | true  | (0.0, 1.0) | "Minimum habitat suitability index."
end


# Rules

# Euler solver rules
@inline applyrule(rule::EulerExponentialGrowth, data, state, args...) = state + state * rule.intrinsicrate * timestep(data)

@inline applyrule(rule::SuitabilityEulerExponentialGrowth, data, state, index, args...) = begin
    intrinsicrate = layer(rule, data, index)
    state + intrinsicrate * state * timestep(data) # dN = rN * dT
    # max(min(state * intrinsicrate, rule.max), rule.min)
end

# TODO: fix and test logistic growth
@inline applyrule(rule::EulerLogisticGrowth, data, state, args...) =
    state + state * rule.intrinsicrate * (oneunit(state) - state / rule.carrycap) * timestep(data) # dN = (1-N/K)rN dT

@inline applyrule(rule::SuitabilityEulerLogisticGrowth, data, state, index, args...) = begin
    intrinsicrate = layer(rule, data, index)
    saturation = intrinsicrate > zero(intrinsicrate) ? (oneunit(state) - state / rule.carrycap) : oneunit(state)
    state + state * saturation * intrinsicrate * timestep(data)
end


# Exact solution rules

@inline applyrule(rule::ExactExponentialGrowth, data, state, args...) =
    state * exp(rule.intrinsicrate * timestep(data))

@inline applyrule(rule::SuitabilityExactExponentialGrowth, data, state, index, args...) = begin
    intrinsicrate = layer(rule, data, index)
    @fastmath state * exp(intrinsicrate * timestep(data))
    # max(min(state * intrinsicrate, rule.max), rule.min)
end

@inline applyrule(rule::ExactLogisticGrowth, data, state, index, args...) = begin
    @fastmath (state * rule.carrycap) / (state + (rule.carrycap - state) * exp(-rule.intrinsicrate * timestep(data)))
end

@inline applyrule(rule::SuitabilityExactLogisticGrowth, data, state, index, args...) = begin
    @inbounds intrinsicrate = layer(rule, data, index)
    # Saturation only applies with positive growth
    if intrinsicrate > zero(intrinsicrate)
        @fastmath (state * rule.carrycap) / (state + (rule.carrycap - state) * exp(-intrinsicrate * timestep(data)))
    else
        @fastmath state * exp(intrinsicrate * timestep(data))
    end
end

@inline applyrule(rule::SuitabilityMask, data, state, index, args...) =
    layer(rule, data, index) >= rule.threshold ? state : zero(state)

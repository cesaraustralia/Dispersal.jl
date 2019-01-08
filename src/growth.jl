# Mixins

@premix @columns struct InstrinsicGrowthRate{GR}
    intrinsicrate::GR = 0.1 | true  | (0.0, 10.0)
end

@premix @columns struct CarryCap{CC}
    carrycap::CC = 100000   | true  | (0.0, 1000000.0)
end

@premix @columns struct Layers{L}
    layers::L    = ()       | false | _
end

# Type declarations

# Euler method solvers

" Simple fixed exponential growth rate solved with Euler method "
@InstrinsicGrowthRate struct EulerExponentialGrowth{} <: AbstractCellModel end

" Simple fixed logistic growth rate solved with Euler method "
@CarryCap @InstrinsicGrowthRate struct EulerLogisticGrowth{} <: AbstractCellModel end

" Exponential growth based on a suitability layer solved with Euler method "
@Layers struct SuitabilityEulerExponentialGrowth{} <: AbstractCellModel end

" Logistic growth based on a suitability layer solved with Euler method "
@CarryCap @Layers struct SuitabilityEulerLogisticGrowth{} <: AbstractCellModel end

# Exact growth solutions

" Simple fixed exponential growth rate using exact solution "
@InstrinsicGrowthRate struct ExactExponentialGrowth{} <: AbstractCellModel end

" Simple fixed logistic growth rate using exact solution "
@CarryCap @InstrinsicGrowthRate struct ExactLogisticGrowth{} <: AbstractCellModel end

" Exponential growth based on a suitability layer using exact solution "
@Layers struct SuitabilityExactExponentialGrowth{} <: AbstractCellModel end

" Logistic growth based on a suitability layer using exact solution "
@CarryCap @Layers struct SuitabilityExactLogisticGrowth{} <: AbstractCellModel end

" Simple suitability layer mask "
@Layers struct SuitabilityMask{ST} <: AbstractCellModel
    # "Minimum habitat suitability index."
    threshold::ST = 0.7 | true  | (0.0, 1.0)
end

# Rules

# Euler solver rules
@inline rule(model::EulerExponentialGrowth, data, state, args...) = state + state * model.intrinsicrate * data.timestep

@inline rule(model::SuitabilityEulerExponentialGrowth, data, state, index, args...) = begin
    state == zero(state) && return state
    intrinsicrate = get_layers(model.layers, index, data.t)
    state + intrinsicrate * state * data.timestep # dN = rN * dT
    # max(min(state * intrinsicrate, model.max), model.min)
end

# TODO: fix and test logistic growth
@inline rule(model::EulerLogisticGrowth, data, state, args...) =
    state + state * model.intrinsicrate * (1 - state / model.carrycap) * data.timestep # dN = (1-N/K)rN dT

@inline rule(model::SuitabilityEulerLogisticGrowth, data, state, index, args...) = begin
    state == zero(state) && return state
    intrinsicrate = get_layers(model.layers, index, data.t)
    saturation =  intrinsicrate > 0 ? (1 - state / model.carrycap) : 1
    state + state * saturation * intrinsicrate * data.timestep
end

# Exact solution rules

@inline rule(model::ExactExponentialGrowth, data, state, args...) = state * exp(model.intrinsicrate * data.timestep)

@inline rule(model::SuitabilityExactExponentialGrowth, data, state, index, args...) = begin
    state == zero(state) && return state
    intrinsicrate = get_layers(model.layers, index, data.t)
    state * exp(intrinsicrate * data.timestep)
    # max(min(state * intrinsicrate, model.max), model.min)
end

@inline rule(model::ExactLogisticGrowth, data, state, args...) =
    (state * model.carrycap) /
    (state + (model.carrycap - state) * exp(-model.intrinsicrate * data.timestep))

@inline rule(model::SuitabilityExactLogisticGrowth, data, state, index, args...) = begin
    state == zero(state) && return state
    intrinsicrate = get_layers(model.layers, index, data.t)
    if intrinsicrate > 0
        (state * model.carrycap) /
        (state + (model.carrycap - state) * exp(-intrinsicrate * data.timestep))
    else
        state * exp(intrinsicrate * data.timestep)
    end
end

@inline rule(model::SuitabilityMask, data, state, index, args...) = begin
    state == zero(state) && return zero(state)
    get_layers(model.layers, index, data.t) >= model.threshold ? state : zero(state)
end

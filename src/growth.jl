# Mixins

@premix @columns struct InstrinsicGrowthRate{GR, TS}
    r::GR = 0.1 | true | (0.0, 10.0)
    timestep::TS = 1 | true | (1, 365.0)
end

@premix @columns struct CarryCap{CC}
    carrycap::CC = 100 | true | (0.0, 10.0)
end

@premix @columns struct Layers{L, TS}
    layers::L = () | false | _
    timestep::TS = 1 | true  | (0.0, 365.0)
end

# Type declarations
# Euler method solvers

" Simple fixed exponential growth rate solved with Euler method "
@InstrinsicGrowthRate struct EulerExponentialGrowth{} <: AbstractModel end

" Simple fixed logistic growth rate solved with Euler method "
@CarryCap @InstrinsicGrowthRate struct EulerLogisticGrowth{} <: AbstractModel end

" Exponential growth based on a suitability layer solved with Euler method "
@Layers struct SuitabilityEulerExponentialGrowth{} <: AbstractModel end

" Logistic growth based on a suitability layer solved with Euler method "
@CarryCap @Layers struct SuitabilityEulerLogisticGrowth{} <: AbstractModel end

# Exact growth solutions

" Simple fixed exponential growth rate using exact solution "
@InstrinsicGrowthRate struct ExactExponentialGrowth{} <: AbstractModel end

" Simple fixed logistic growth rate using exact solution "
@CarryCap @InstrinsicGrowthRate struct ExactLogisticGrowth{} <: AbstractModel end

" Exponential growth based on a suitability layer using exact solution "
@Layers struct SuitabilityExactExponentialGrowth{} <: AbstractModel end

" Logistic growth based on a suitability layer using exact solution "
@CarryCap @Layers struct SuitabilityExactLogisticGrowth{} <: AbstractModel end

" Simple suitability layer mask "
@Layers struct SuitabilityMask{ST} <: AbstractModel
    # "Minimum habitat suitability index."
    threshold::ST = 0.7 | true  | (0.0, 1.0)
end

# Rules
# euler solver rules
@inline rule(model::EulerExponentialGrowth, data, state, args...) = state + state * model.r

@inline rule(model::SuitabilityEulerExponentialGrowth, data, state, index, args...) = begin
    state == zero(state) && return state
    r = get_layers(model.layers, index, data.t)
    state + r * state * model.timestep # dN = rN * dT
    # max(min(state * r, model.max), model.min)
end

# TODO: fix and test logistic growth
@inline rule(model::EulerLogisticGrowth, data, state, args...) =
    state + state * model.r * (1 - state / model.carrycap) * model.timestep # dN = (1-N/K)rN dT

@inline rule(model::SuitabilityEulerLogisticGrowth, data, state, index, args...) = begin
    state == zero(state) && return state
    r = get_layers(model.layers, index, data.t)
    saturation =  r > 0 ? (1 - state / model.carrycap) : 1
    state + state * saturation * r * model.timestep
end

# exact solution rules

@inline rule(model::ExactExponentialGrowth, data, state, args...) = state * exp(model.r * model.timestep)

@inline rule(model::SuitabilityExactExponentialGrowth, data, state, index, args...) = begin
    state == zero(state) && return state
    r = get_layers(model.layers, index, data.t)
    state * exp(r * model.timestep)
    # max(min(state * r, model.max), model.min)
end

@inline rule(model::ExactLogisticGrowth, data, state, args...) =
    (state * model.carrycap) /
    (state + (model.carrycap - state) * exp(-model.r * model.timestep))

@inline rule(model::SuitabilityExactLogisticGrowth, data, state, index, args...) = begin
    state == zero(state) && return state
    r = get_layers(model.layers, index, data.t)
    if r > 0
        (state * model.carrycap) /
        (state + (model.carrycap - state) * exp(-r * model.timestep))
    else
        state * exp(r * model.timestep)
    end
end

@inline rule(model::SuitabilityMask, data, state, index, args...) = begin
    state == zero(state) && return zero(state)
    get_layers(model.layers, index, data.t) >= model.threshold ? state : zero(state)
end

# Mixins

@premix @columns struct InstrinsicGrowthRate{GR, TS}
    intinsicrate::GR = 0.1 | true | (0.0, 10.0)
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

" Simple fixed exponential growth intinsicrate solved with Euler method "
@InstrinsicGrowthRate struct EulerExponentialGrowth{} <: AbstractModel end

" Simple fixed logistic growth intinsicrate solved with Euler method "
@CarryCap @InstrinsicGrowthRate struct EulerLogisticGrowth{} <: AbstractModel end

" Exponential growth based on a suitability layer solved with Euler method "
@Layers struct SuitabilityEulerExponentialGrowth{} <: AbstractModel end

" Logistic growth based on a suitability layer solved with Euler method "
@CarryCap @Layers struct SuitabilityEulerLogisticGrowth{} <: AbstractModel end

# Exact growth solutions

" Simple fixed exponential growth intinsicrate using exact solution "
@InstrinsicGrowthRate struct ExactExponentialGrowth{} <: AbstractModel end

" Simple fixed logistic growth intinsicrate using exact solution "
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
@inline rule(model::EulerExponentialGrowth, data, state, args...) = state + state * model.intinsicrate * model.timestep

@inline rule(model::SuitabilityEulerExponentialGrowth, data, state, index, args...) = begin
    state == zero(state) && return state
    intinsicrate = get_layers(model.layers, index, data.t)
    state + intinsicrate * state * model.timestep # dN = rN * dT
    # max(min(state * intinsicrate, model.max), model.min)
end

# TODO: fix and test logistic growth
@inline rule(model::EulerLogisticGrowth, data, state, args...) =
    state + state * model.intinsicrate * (1 - state / model.carrycap) * model.timestep # dN = (1-N/K)rN dT

@inline rule(model::SuitabilityEulerLogisticGrowth, data, state, index, args...) = begin
    state == zero(state) && return state
    intinsicrate = get_layers(model.layers, index, data.t)
    saturation =  intinsicrate > 0 ? (1 - state / model.carrycap) : 1
    state + state * saturation * intinsicrate * model.timestep
end

# exact solution rules

@inline rule(model::ExactExponentialGrowth, data, state, args...) = state * exp(model.intinsicrate * model.timestep)

@inline rule(model::SuitabilityExactExponentialGrowth, data, state, index, args...) = begin
    state == zero(state) && return state
    intinsicrate = get_layers(model.layers, index, data.t)
    state * exp(intinsicrate * model.timestep)
    # max(min(state * intinsicrate, model.max), model.min)
end

@inline rule(model::ExactLogisticGrowth, data, state, args...) =
    (state * model.carrycap) /
    (state + (model.carrycap - state) * exp(-model.intinsicrate * model.timestep))

@inline rule(model::SuitabilityExactLogisticGrowth, data, state, index, args...) = begin
    state == zero(state) && return state
    intinsicrate = get_layers(model.layers, index, data.t)
    if intinsicrate > 0
        (state * model.carrycap) /
        (state + (model.carrycap - state) * exp(-intinsicrate * model.timestep))
    else
        state * exp(intinsicrate * model.timestep)
    end
end

@inline rule(model::SuitabilityMask, data, state, index, args...) = begin
    state == zero(state) && return zero(state)
    get_layers(model.layers, index, data.t) >= model.threshold ? state : zero(state)
end

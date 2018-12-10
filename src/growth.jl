# Mixins

@premix @columns struct GrowthRate{GR}
    growthrate::GR = 0.1 | true | (0.0, 10.0)
end

@premix @columns struct CarryCap{CC}
    carrycap::CC = 100 | true | (0.0, 10.0)
end

@premix @columns struct Layers{L}
    layers::L = () | false | _
end

# Type declarations

" Simple fixed exponential growth rate "
@GrowthRate struct ExponentialGrowth{} <: AbstractModel end

" Simple fixed logistic growth rate "
@CarryCap @GrowthRate struct LogisticGrowth{} <: AbstractModel end

" Exponential growth based on a suitability layer "
@Layers struct SuitabilityExponentialGrowth{} <: AbstractModel end

" Logistic growth based on a suitability layer "
@CarryCap @Layers struct SuitabilityLogisticGrowth{} <: AbstractModel end

" Simple suitability layer mask "
@Layers struct SuitabilityMask{ST} <: AbstractModel
    # "Minimum habitat suitability index."
    threshold::ST = 0.7 | true  | (0.0, 1.0)
end

# Rules

@inline rule(model::ExponentialGrowth, data, state, args...) = state + state * model.growthrate

@inline rule(model::SuitabilityMask, data, state, index, args...) = begin
    state == zero(state) && return zero(state)
    get_layers(model.layers, index, data.t) >= model.threshold ? state : zero(state)
end

@inline rule(model::SuitabilityExponentialGrowth, data, state, index, args...) = begin
    state == zero(state) && return state
    growthrate = get_layers(model.layers, index, data.t)
    state + growthrate * state
    # max(min(state * growthrate, model.max), model.min)
end

# TODO: fix and test logistic growth
@inline rule(model::LogisticGrowth, data, state, args...) =
    state + state * model.growthrate * (1 - state / model.carrycap)

@inline rule(model::SuitabilityLogisticGrowth, data, state, index, args...) = begin
    state == zero(state) && return state
    growthrate = get_layers(model.layers, index, data.t)
    saturation =  growthrate > 0 ? (1 - state / model.carrycap) : 1
    state + state * saturation * growthrate
end

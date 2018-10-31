@mix @columns struct GrowthRate{GR}
    growthrate::GR = 1.1 | true | (0.0, 10.0)
end

@mix @columns struct MinMax{M}
    min::M = 0.0       | false | _
    max::M = 1000000.0 | false | _
end

@mix @columns struct CarryCap{CC}
    carrycap::CC = 1.1 | true | (0.0, 10.0)
end


" Simple fixed exponential growth rate "
@GrowthRate struct ExponentialGrowth{} <: AbstractModel end

" Simple fixed logistic growth rate "
@GrowthRate @CarryCap struct LogisticGrowth{} <: AbstractModel end


" Exponential growth based on a suitability layer "
@MinMax struct SuitabilityExponentialGrowth{} <: AbstractModel end

" Logistic growth based on a suitability layer "
@MinMax @CarryCap struct SuitabilityLogisticGrowth{} <: AbstractModel end

SuitabilityExponentialGrowth(init::ScalableMatrix, args...) = 
    SuitabilityExponentialGrowth(init.min, init.max, args...)

SuitabilityLogisticGrowth(init::ScalableMatrix, args...) = 
    SuitabilityLogisticGrowth(init.min, init.max, args...)

" Simple suitability layer mask "
@columns struct SuitabilityMask{ST} <: AbstractModel 
    # "Minimum habitat suitability index."
    suitability_threshold::ST = 0.7 | true  | (0.0, 1.0)
end


rule(model::ExponentialGrowth, state, args...) = state * model.growthrate
rule(model::LogisticGrowth, state, args...) = state + state * model.rate * (1 - state / model.carrycap) 

rule(model::SuitabilityMask, state, index, t, source, dest, layers, args...) = begin
    state == zero(state) && return state
    suit = suitability(layers, index, t)
    return suit >= model.suitability_threshold ? state : zero(state)
end

rule(model::SuitabilityExponentialGrowth, state, index, t, source, dest, layers, args...) = begin
    state == zero(state) && return state
    growthrate = suitability(layers, index, t)
    max(min(state * growthrate, model.max), model.min)
end

rule(model::SuitabilityLogisticGrowth, state, index, t, source, dest, layers, args...) = begin
    state == zero(state) && return state
    growthrate = suitability(layers, index, t)
    state1 = state + state * growthrate * (1 - state / model.carrycap) 
    max(min(state1 , model.max), model.min)
end


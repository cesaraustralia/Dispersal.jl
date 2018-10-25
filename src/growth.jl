@mix struct GrowthRate 
    growthrate::R = 1.1 | true | (0.0, 10.0)
end

@mix struct MinMax
    min::M       = 0.0       | false | _
    max::M       = 1000000.0 | false | _
end

@mix struct CarryCap
    carrycap::R = 1.1 | true | (0.0, 10.0)
end


" Simple fixed exponential growth rate "
@GrowthRate struct ExponentialGrowth{R} <: AbstractModel end

" Simple fixed logistic growth rate "
@GrowthRate @CarryCap struct LogisticGrowth{R} <: AbstractModel end


" Exponential growth based on a suitability layer "
@MinMax struct SuitabilityExponentialGrowth{M,TS} <: AbstractModel end

" Logistic growth based on a suitability layer "
@MinMax @CarryCap struct SuitabilityLogisticGrowth{M,TS} <: AbstractModel end

SuitabilityExponentialGrowth(init::ScalableMatrix, args...) = 
    SuitabilityExponentialGrowth(init.min, init.max, args...)

SuitabilityLogisticGrowth(init::ScalableMatrix, args...) = 
    SuitabilityLogisticGrowth(init.min, init.max, args...)

" Simple suitability layer mask "
@columns struct SuitabilityMask{ST} <: AbstractModel 
    # "Minimum habitat suitability index."
    suitability_threshold::ST = 0.1 | true  | (0.0, 1.0)
end


rule(model::ExponentialGrowth, state, args...) = state * model.growthrate
rule(model::LogisticGrowth, state, args...) = state * model.growthrate

rule(model::SuitabilityMask, state, row, col, t, source, dest, layers, args...) = begin
    state == zero(state) && return state
    suit = suitability(layers, (row, col), t)
    return suit >= model.suitability_threshold ? state : zero(state)
end

rule(model::SuitabilityExponentialGrowth, state, row, col, t, source, dest, layers, args...) = begin
    state == zero(state) && return state
    growthrate = suitability(layers, (row, col), t)
    state1 = state + state * growthrate * (1 - state / model.carrycap) 
    max(min(state1 , model.max), model.min)
end

rule(model::SuitabilityLogisticGrowth, state, row, col, t, source, dest, layers, args...) = begin
    state == zero(state) && return state
    growthrate = suitability(layers, (row, col), t)
    max(min(state * growthrate, model.max), model.min)
end

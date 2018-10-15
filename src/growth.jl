" Simple linear growth rate "
@columns struct FixedGrowth{R} <: AbstractModel 
    growthrate::R = 1.1 | true | (0.0, 10.0)
end

" Growth model based on suitability layer "
@columns struct SuitabilityGrowth{M,TS} <: AbstractModel 
    min::M       = 0.0       | false | _
    max::M       = 1000000.0 | false | _
    timestep::TS = 1.0       | true  | (0.0, 10.0) 
end

SuitabilityGrowth(init::ScalableMatrix, timestep=1/30.0) = 
    SuitabilityGrowth(init.min, init.max, timestep)

" Growth model based on suitability layer "
@columns struct SuitabilityMask{ST} <: AbstractModel 
    # "Minimum habitat suitability index."
    suitability_threshold::ST = 0.1 | true  | (0.0, 1.0)
end

" Grow the population at a fixed rate "
rule(model::FixedGrowth, state, args...) = state * model.growthrate

rule(model::SuitabilityMask, state, row, col, t, source, dest, layers, args...) = begin
    state == zero(state) && return state
    suit = suitability(layers, (row, col), t)
    return suit >= model.suitability_threshold ? state : zero(state)
end

rule(model::SuitabilityGrowth, state, row, col, t, source, dest, layers, args...) = begin
    state == zero(state) && return state
    growthrate = exp(suitability(layers, (row, col), t) * model.timestep)
    max(min(state * growthrate, model.max), model.min)
end

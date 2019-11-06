# Mixins

@premix @columns struct Timestep{TS,S}
    # Field            | Default | Flatn | Bounds      | Description
    timestep::TS       | nothing | false | _           | "Timestep converted from sim data. Needs to be separate from rate for DateTime"
    nsteps::S          | 1.0     | false | _           | "The exact nsteps timestep, updated by precalcrule"
end

@premix @Timestep struct InstrinsicGrowthRate{GR}
    intrinsicrate::GR  | 0.1     | true  | (0.0, 10.0) | "Intrinsic rate of growth per timestep"
end

@premix @Timestep struct Layers{L,TI}
    layer::L          | nothing  | false | _           | "Data layer"
    timeinterp::TI    | 1        | false | _           | "Precalculated interpolation indices"
end

@premix @columns struct CarryCap{CC}
    carrycap::CC      | 100000.0 | true  | (0.0, 1e9)  | "Carrying capacity for each cell. Not currently scaled by area."
end


"""
Extends CellRule for rules of growth dynamics

For best performance these should be chained with other
CellRule or following an NeighborhoodRule.
"""
abstract type GrowthRule <: CellRule end


"""
Extends GrowthRule for growth rules using a heterogenous
growth rate layer. 

[GrowthMaps.jl](http://github.com/cesaraustralia/GrowthMaps.jl)
can produce these growth maps from environmental data.
"""
abstract type GrowthMapRule <: GrowthRule end

layer(rule::GrowthMapRule) = rule.layer 
timeinterp(rule::GrowthMapRule) = rule.timeinterp

DynamicGrids.precalcrules(rule::GrowthMapRule, data) = begin
    if :timestep in fieldnames(typeof(rule))
        rule = precaltimestep(rule, data)
    end
    DynamicGrids.precalcrules(layer(rule), rule, data)
end
DynamicGrids.precalcrules(::AbstractMatrix, rule::GrowthMapRule, data) = rule
DynamicGrids.precalcrules(::AbstractArray{<:Any,3}, rule::GrowthMapRule, data) =
    @set rule.timeinterp = precalc_time_interpolation(layer(rule), rule, data)
DynamicGrids.precalcrules(rule::GrowthRule, data) = 
    if :timestep in fieldnames(typeof(rule))
        precaltimestep(rule, data)
    else
        rule
    end

precaltimestep(rule, data) = precaltimestep(rule.timestep, rule, data)
precaltimestep(ruletimestep::DatePeriod, rule, data) = 
    @set rule.nsteps = currenttimestep(data) / Millisecond(ruletimestep) 
precaltimestep(ruletimestep::Nothing, rule, data) = @set rule.nsteps = 1 
precaltimestep(ruletimestep, rule, data) = 
    @set rule.nsteps = timestep(data) / ruletimestep


# Exact growth solutions

"""
Simple fixed exponential growth rate using exact solution.
$(FIELDDOCTABLE)
"""
@InstrinsicGrowthRate struct ExactExponentialGrowth{} <: GrowthRule end

@inline applyrule(rule::ExactExponentialGrowth, data, state, args...) =
    state * exp(rule.intrinsicrate * rule.nsteps)

"""
Simple fixed logistic growth rate using exact solution
$(FIELDDOCTABLE)
"""
@CarryCap @InstrinsicGrowthRate struct ExactLogisticGrowth{} <: GrowthRule end

@inline applyrule(rule::ExactLogisticGrowth, data, state, index, args...) = begin
    @fastmath (state * rule.carrycap) / (state + (rule.carrycap - state) * 
                                         exp(-rule.intrinsicrate * rule.nsteps))
end

"""
Exponential growth based on a growth rate layer using exact solution.
$(FIELDDOCTABLE)
"""
@Layers struct ExactExponentialGrowthMap{} <: GrowthMapRule end

@inline applyrule(rule::ExactExponentialGrowthMap, data, state, index, args...) = begin
    intrinsicrate = layer(rule, data, index)
    @fastmath state * exp(intrinsicrate * rule.nsteps)
end

"""
Logistic growth based on a growth rate layer, using exact solution.

Saturation only applies with positive growth
$(FIELDDOCTABLE)
"""
@CarryCap @Layers struct ExactLogisticGrowthMap{} <: GrowthMapRule end

@inline applyrule(rule::ExactLogisticGrowthMap, data, state, index, args...) = begin
    @inbounds intrinsicrate = layer(rule, data, index)
    if intrinsicrate > zero(intrinsicrate)
        @fastmath (state * rule.carrycap) / (state + (rule.carrycap - state) * 
                                             exp(-intrinsicrate * rule.nsteps))
    else
        @fastmath state * exp(intrinsicrate * rule.nsteps)
    end
end

"""
Simple layer mask. Values below a certain threshold are replaced with zero.
$(FIELDDOCTABLE)
"""
@Layers struct MaskGrowthMap{ST} <: GrowthMapRule
    threshold::ST | 0.5 | true  | (0.0, 1.0) | "Minimum viability index."
end

@inline applyrule(rule::MaskGrowthMap, data, state, index, args...) =
    layer(rule, data, index) >= rule.threshold ? state : zero(state)

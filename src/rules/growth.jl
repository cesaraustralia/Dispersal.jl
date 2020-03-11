# Mixins

@mix @columns struct CarryCap{CC}
    # Field           | Default  | Flat  | Bounds      | Description
    carrycap::CC      | 100000.0 | true  | (0.0, 1e9)  | "Carrying capacity for each cell. Not currently scaled by area."
end

@mix @columns struct Timestep{TS,NS}
    timestep::TS      | nothing  | false | _           | "Timestep converted from sim data. Needs to be separate from rate for DateTime"
    nsteps::NS        | 1.0      | false | _           | "The exact nsteps timestep, updated by precalcrule"
end

@mix @Timestep struct InstrinsicGrowthRate{GR}
    intrinsicrate::GR | 0.1      | true  | (0.0, 10.0) | "Intrinsic rate of growth per timestep"
end

precalctimestep(rule, data) = precalctimestep(rule.timestep, rule, data)
precalctimestep(ruletimestep::DatePeriod, rule, data) =
    @set rule.nsteps = typeof(rule.nsteps)(currenttimestep(data) / Millisecond(ruletimestep))
precalctimestep(ruletimestep::Nothing, rule, data) = @set rule.nsteps = oneunit(rule.nsteps)
precalctimestep(ruletimestep, rule, data) =
    @set rule.nsteps = typeof(rule.nsteps)(timestep(data) / ruletimestep)

"""
Extends CellRule for rules of growth dynamics

For best performance these should be chained with other
CellRule or following an NeighborhoodRule.
"""
abstract type GrowthRule{R,W} <: CellRule{R,W} end


"""
Extends GrowthRule for growth rules using a heterogenous
growth rate layer.

[GrowthMaps.jl](http://github.com/cesaraustralia/GrowthMaps.jl)
can produce these growth maps from environmental data.
"""
abstract type GrowthMapRule{R,W} <: GrowthRule{R,W} end

layer(rule::GrowthMapRule) = rule.layer
timeinterp(rule::GrowthMapRule) = rule.timeinterp

DynamicGrids.precalcrules(rule::GrowthMapRule, data) = begin
    if :timestep in fieldnames(typeof(rule))
        rule = precalctimestep(rule, data)
    end
    precalclayer(layer(rule), rule, data)
end

DynamicGrids.precalcrules(rule::GrowthRule, data) =
    if :timestep in fieldnames(typeof(rule))
        precalctimestep(rule, data)
    else
        rule
    end


# Exact growth solutions

"""
Simple fixed exponential growth rate using exact solution.
$(FIELDDOCTABLE)
"""
@InstrinsicGrowthRate struct ExactExponentialGrowth{R,W} <: GrowthRule{R,W} end

@inline applyrule(rule::ExactExponentialGrowth, data, state, args...) =
    state * exp(rule.intrinsicrate * rule.nsteps)

"""
Simple fixed logistic growth rate using exact solution
$(FIELDDOCTABLE)
"""
@InstrinsicGrowthRate @CarryCap struct ExactLogisticGrowth{R,W} <: GrowthRule{R,W} end

@inline applyrule(rule::ExactLogisticGrowth, data, state, index, args...) = begin
    @fastmath (state * rule.carrycap) / (state + (rule.carrycap - state) *
                                         exp(-rule.intrinsicrate * rule.nsteps))
end

"""
Exponential growth based on a growth rate layer using exact solution.
$(FIELDDOCTABLE)
"""
@Timestep @Layers struct ExactExponentialGrowthMap{R,W} <: GrowthMapRule{R,W} end

@inline applyrule(rule::ExactExponentialGrowthMap, data, state, index, args...) = begin
    intrinsicrate = layer(rule, data, index)
    @fastmath state * exp(intrinsicrate * rule.nsteps)
end

"""
Logistic growth based on a growth rate layer, using exact solution.

Saturation only applies with positive growth
$(FIELDDOCTABLE)
"""
@Layers @Timestep @CarryCap struct ExactLogisticGrowthMap{R,W} <: GrowthMapRule{R,W} end

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
@Layers @Timestep struct MaskGrowthMap{R,W,ST} <: GrowthMapRule{R,W}
    threshold::ST | 0.5 | true  | (0.0, 1.0) | "Minimum viability index."
end

@inline applyrule(rule::MaskGrowthMap, data, state, index, args...) =
    layer(rule, data, index) >= rule.threshold ? state : zero(state)

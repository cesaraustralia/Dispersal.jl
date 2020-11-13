
const CARRYCAP      = Param(100000.0; bounds=(0.0, 10.0))
const INTRINSICRATE = Param(0.1,      bounds=(0.0, 10.0))
const THRESHOLD     = Param(0.5;      bounds=(0.0, 1.0))

"""
Extends CellRule for rules of growth dynamics

For best performance these should be chained with other
CellRule or following an NeighborhoodRule.
"""
abstract type GrowthRule{R,W} <: CellRule{R,W} end


"""
Extends GrowthRule for growth rules using heterogenous
growth rate data.
a
[GrowthMaps.jl](http://github.com/cesaraustralia/GrowthMaps.jl)
can produce these growth maps from environmental data.

$(FIELDDOCTABLE)
"""
abstract type GrowthMapRule{R,W} <: GrowthRule{R,W} end


"""
Simple fixed exponential growth rate using exact solution.

$(FIELDDOCTABLE)
"""
struct ExactExponentialGrowth{R,W,GR,TS,S} <: GrowthRule{R,W}
    "Intrinsic rate of growth per timestep"
    intrinsicrate::GR
    "Timestep used in formulation"
    timestep::TS
    "The fractional number of rule timesteps in the current simulation timestep"
    nsteps::S
end
ExactExponentialGrowth{R,W}(;
    intrinsicrate=INTRINSICRATE,
    timestep=nothing,
    nsteps=1.0,
) where {R,W} = ExactExponentialGrowth{R,W}(intrinsicrate, timestep, nsteps)

# DynamicGrids.jl interface

precalcrule(rule::ExactExponentialGrowth, data) = precalc_timestep(rule, data)

@inline function applyrule(data, rule::ExactExponentialGrowth, population, cellindex)
    population > zero(population) || return zero(population)
    population * exp(rule.intrinsicrate * rule.nsteps)
end


"""
Simple fixed logistic growth rate using exact solution

$(FIELDDOCTABLE)
"""
struct ExactLogisticGrowth{R,W,CC,GR,TS,S} <: GrowthRule{R,W}
    "Carrying capacity for each cell. Not currently scaled by area."
    carrycap::CC
    "Intrinsic rate of growth per timestep"
    intrinsicrate::GR
    "Timestep used in formulation"
    timestep::TS
    "The fractional number of rule timesteps in the current simulation timestep"
    nsteps::S
end
ExactLogisticGrowth{R,W}(;
    carrycap=CARRYCAP,
    intrinsicrate=INTRINSICRATE,
    timestep=nothing,
    nsteps=1.0,
    ) where {R,W} = ExactLogisticGrowth{R,W}(carrycap, intrinsicrate, timestep, nsteps)

# DynamicGrids.jl interface

precalcrule(rule::ExactLogisticGrowth, data) = precalc_timestep(rule, data)

@inline function applyrule(data, rule::ExactLogisticGrowth, population, cellindex)
    population > zero(population) || return zero(population)
    @fastmath (population * rule.carrycap) / (population + (rule.carrycap - population) *
                                         exp(-rule.intrinsicrate * rule.nsteps))
end



"""
Exponential growth based on a growth rate data, using exact solution.

$(FIELDDOCTABLE)
"""
struct ExactExponentialGrowthMap{R,W,AK<:Val,AT,TS,S} <: GrowthMapRule{R,W}
    "Key for aux data"
    auxkey::AK
    "Precalculated time interpolation index for aux data"
    auxtimeindex::AT
    "Timestep used in formulation"
    timestep::TS
    "The fractional number of rule timesteps in the current simulation timestep"
    nsteps::S
end
ExactExponentialGrowthMap{R,W}(;
    auxkey,
    auxtimeindex=1,
    timestep=nothing,
    nsteps=1.0,
) where {R,W} = ExactExponentialGrowthMap{R,W}(auxkey, auxtimeindex, timestep, nsteps)

# DynamicGrids.jl interface

function precalcrule(rule::ExactExponentialGrowthMap, data)
    rule = precalc_timestep(rule, data)
    precalc_auxtimeindex(aux(data, rule.auxkey), rule, data)
end

@inline function applyrule(data, rule::ExactExponentialGrowthMap, population, cellindex)
    population > zero(population) || return zero(population)
    intrinsicrate = auxval(data, rule.auxkey, cellindex..., rule.auxtimeindex) 
    @fastmath population * exp(intrinsicrate * rule.nsteps)
end


"""
Logistic growth based on a growth rate layer, using exact solution.

Saturation only applies with positive growth
$(FIELDDOCTABLE)
"""
struct ExactLogisticGrowthMap{R,W,AK<:Val,AT,CC,TS,S} <: GrowthMapRule{R,W}
    "Key for aux layer"
    auxkey::AK
    "Precalculated time interpolation index for aux data"
    auxtimeindex::AT
    "Carrying capacity for each cell. Not currently scaled by area."
    carrycap::CC
    "Timestep used in formulation"
    timestep::TS
    "The fractional number of rule timesteps in the current simulation timestep"
    nsteps::S
end
ExactLogisticGrowthMap{R,W}(;
    auxkey,
    auxtimeindex=:1,
    carrycap=CARRYCAP,
    timestep=nothing,
    nsteps=1.0,
) where {R,W} = ExactLogisticGrowthMap{R,W}(auxkey, auxtimeindex, carrycap, timestep, nsteps)

# DynamicGrids.jl interface

function precalcrule(rule::ExactLogisticGrowthMap, data)
    rule = precalc_timestep(rule, data)
    precalc_auxtimeindex(aux(data, rule.auxkey), rule, data)
end

@inline function applyrule(data, rule::ExactLogisticGrowthMap, population, cellindex)
    population > zero(population) || return zero(population)
    intrinsicrate = auxval(data, rule.auxkey, cellindex..., rule.auxtimeindex) 
    if intrinsicrate > zero(intrinsicrate)
        @fastmath (population * rule.carrycap) / (population + (rule.carrycap - population) *
                                             exp(-intrinsicrate * rule.nsteps))
    else
        @fastmath population * exp(intrinsicrate * rule.nsteps)
    end
end


"""
Simple layer mask. Values below a certain threshold are replaced with zero.

$(FIELDDOCTABLE)
"""
struct MaskGrowthMap{R,W,AK<:Val,AT,Th} <: GrowthMapRule{R,W}
    "Key for aux data"
    auxkey::AK
    "Precalculated time interpolation index for aux data"
    auxtimeindex::AT
    "Minimum viability threshold."
    threshold::Th
end
MaskGrowthMap{R,W}(;
    auxkey,
    auxtimeindex=1,
    threshold=THRESHOLD,
) where {R,W} = MaskGrowthMap{R,W}(auxkey, auxtimeindex, threshold)

# DynamicGrids.jl interface

function precalcrule(rule::MaskGrowthMap, data)
    precalc_auxtimeindex(aux(data, rule.auxkey), rule, data)
end

@inline function applyrule(data, rule::MaskGrowthMap, population, cellindex)
    intrinsicrate = auxval(data, rule.auxkey, cellindex..., rule.auxtimeindex) 
    intrinsicrate >= rule.threshold ? population : zero(population)
end

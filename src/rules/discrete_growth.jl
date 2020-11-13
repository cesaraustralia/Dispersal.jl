const INTRINSICRATE_PARAM = Param(0.1, bounds=(0.0, 10.0))

"""
Simple discrete growth rate.
"""
struct DiscreteGrowth{R,W,GR} <: GrowthRule{R,W}
    "Intrinsic rate of growth per timestep"
    intrinsicrate::GR
end
function DiscreteGrowth{R,W}(; intrinsicrate=INTRINSICRATE_PARAM) where {R,W}
    DiscreteGrowth{R,W}(intrinsicrate)
end

@inline function applyrule(data, rule::DiscreteGrowth, population, args...)
    population > zero(population) || return zero(population)
    rule.intrinsicrate * population
end

struct CappedDiscreteGrowth{R,W,GR,CC} <: GrowthRule{R,W}
    "Intrinsic rate of growth per timestep"
    intrinsicrate::GR
    "Carrying capacity"
    carrycap::CC
end
function CappedDiscreteGrowth{R,W}(; 
    intrinsicrate=INTRINSICRATE_PARAM,
    carrycap=CARRYCAP_PARAM,
) where {R,W}
    CappedDiscreteGrowth{R,W}(intrinsicrate, carrycap)
end

@inline function applyrule(data, rule::CappedDiscreteGrowth, pops, args...)
    carrycap_effect = (1 - sum(pops) / rule.carrycap)
    map(values(pops),  rule.intrinsicrate) do pop, ir
        ir * pop * carrycap_effect
    end
end

"""
Simple discrete growth rate.
"""
struct DiscreteGrowthMap{R,W,RK,AT} <: GrowthMapRule{R,W}
    "key for carrycap layer"
    ratekey::RK
    "Precalculated time interpolation index for aux data"
    auxtimeindex::AT
end
function DiscreteGrowthMap{R,W}(; ratekey, auxtimeindex=1) where {R,W}
    DiscreteGrowthMap{R,W}(ratekey, auxtimeindex)
end

@inline function applyrule(data, rule::DiscreteGrowthMap, population, index, args...)
    population > zero(population) || return zero(population)
    intrinsicrate = auxval(data, rule.ratekey, index..., rule.auxtimeindex) 
    @fastmath intrinsicrate * population
end

struct LogisticGrowthMap3{R,W,RK,CK,AT,TS,S} <: GrowthMapRule{R,W}
    "Key for growth rate layer"
    ratekey::RK
    "key for carrycap layer"
    carrycapkey::CK
    "Precalculated time interpolation index for aux data"
    auxtimeindex::AT
    "Timestep used in formulation"
    timestep::TS
    "The fractional number of rule timesteps in the current simulation timestep"
    nsteps::S
end
function LogisticGrowthMap3{R,W}(;
    ratekey,
    carrycapkey,
    auxtimeindex=1,
    timestep=nothing,
    nsteps=1
) where {R,W}
    LogisticGrowthMap3{R,W}(ratekey, carrycapkey, auxtimeindex, timestep, nsteps)
end

@inline function applyrule(data, rule::LogisticGrowthMap3, population, index, args...)
    population > zero(population) || return zero(population)
    # I'll think of a cleaner way to do this part
    intrinsicrate = aux(data)[unwrap(rule.ratekey)][index..., rule.auxtimeindex]
    carrycap = aux(data)[unwrap(rule.carrycapkey)][index..., rule.auxtimeindex]

    if intrinsicrate > zero(intrinsicrate)
        @fastmath (population * carrycap) / (population + (carrycap - population) *
                                             exp(-intrinsicrate * rule.nsteps))
    else
        @fastmath population * exp(intrinsicrate * rule.nsteps)
    end
end

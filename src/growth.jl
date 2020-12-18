
const CARRYCAP_PARAM      = Param(100000.0; bounds=(0.0, 10.0))
const INTRINSICRATE_PARAM = Param(0.1,      bounds=(0.0, 10.0))
const THRESHOLD_PARAM     = Param(0.5;      bounds=(0.0, 1.0))

"""
Extends CellRule for rules of growth dynamics

For best performance these should be chained with other
CellRule or following an NeighborhoodRule.
"""
abstract type GrowthRule{R,W} <: CellRule{R,W} end

"""
Exponential growth based on a growth rate data, using exact solution.
"""
struct ExponentialGrowth{R,W,GR,TS,S} <: GrowthRule{R,W}
    "Key for aux data or single rate"
    rate::GR
    "Precalculated time interpolation index for aux data"
    timestep::TS
    "The fractional number of rule timesteps in the current simulation timestep"
    nsteps::S
end
function ExponentialGrowth{R,W}(; 
    rate=INTRINSICRATE_PARAM, 
    timestep=nothing, 
    nsteps=1.0
) where {R,W}
    ExponentialGrowth{R,W}(rate, timestep, nsteps)
end

precalcrule(rule::ExponentialGrowth, data) = precalc_timestep(rule, data)

@inline function applyrule(data, rule::ExponentialGrowth, population, I)
    population > zero(population) || return zero(population)
    intrinsicrate = get(data, rule.rate, I...)
    @fastmath population * exp(intrinsicrate * rule.nsteps)
end

"""
Logistic growth based on a growth rate layer, using exact solution.

Saturation only applies with positive growth
"""
struct LogisticGrowth{R,W,GR,CC,TS,S} <: GrowthRule{R,W}
    "Key for aux data or single rate"
    rate::GR
    "Carrying capacity for each cell"
    carrycap::CC
    "Timestep used in the formulation"
    timestep::TS
    "The fractional number of rule timesteps in the current simulation timestep"
    nsteps::S
end
function LogisticGrowth{R,W}(;
    rate=INTRINSICRATE_PARAM, 
    carrycap=CARRYCAP_PARAM, 
    timestep=nothing, nsteps=1.0,
) where {R,W}
    LogisticGrowth{R,W}(rate, carrycap, timestep, nsteps)
end

precalcrule(rule::LogisticGrowth, data) = precalc_timestep(rule, data)

@inline function applyrule(data, rule::LogisticGrowth, population, I)
    population > zero(population) || return zero(population)
    intrinsicrate = get(data, rule.rate, I...) 
    carrycap = get(data, rule.carrycap, I...) 

    if intrinsicrate > zero(intrinsicrate)
        @fastmath (population * carrycap) / (population + (carrycap - population) *
                  exp(-intrinsicrate * rule.nsteps))
    else
        @fastmath population * exp(intrinsicrate * rule.nsteps)
    end
end

"""
Simple threshold mask. Values below a certain threshold are replaced with zero.
"""
struct ThresholdGrowth{R,W,GR,Th} <: GrowthRule{R,W}
    "Key for aux data or single rate"
    rate::GR
    "Minimum viability threshold."
    threshold::Th
end
function ThresholdGrowth{R,W}(; 
    rate=INTRINSICRATE_PARAM, 
    threshold=THRESHOLD_PARAM,
) where {R,W}
    ThresholdGrowth{R,W}(rate, threshold)
end

@inline function applyrule(data, rule::ThresholdGrowth, population, I)
    intrinsicrate = get(data, rule.rate, I...)
    threshold = get(data, rule.threshold, I...)
    intrinsicrate >= threshold ? population : zero(population)
end

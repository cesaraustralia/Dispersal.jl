const LC50_PARAM = Param(10.0; bounds=(0.0, 1e9))
const HILLCOEFFICIENT_PARAM = Param(0.1; bounds=(0.0, 1e9))
const EXPOSURE_PARAM = Param(100.0, bounds=(0.0, 1e9))

"""
Extends CellRule for rules of survival effect

For best performance these should be chained with other
CellRule or following an NeighborhoodRule.
"""
abstract type SurvRule{R,W} <: CellRule{R,W} end

"""
Extends SurvRule for survival rules using a heterogenous
concentration layer.
"""

abstract type SurvMapRule{R,W} <: SurvRule{R,W} end

function DynamicGrids.precalcrules(rule::SurvMapRule, data)
    precalc_auxtimeindex(aux(data, rule.exposurekey), rule, data)
end

# Exact survival solutions

"""
Simple fixed survival loglogistic using exact solution.

Cumulative function of loglogistic: F(x; α, β) =  1/(1 + (x/ α)^{-β} ),
    where α>0 is the scale and β>0 is the shape.

Dose-response log-logistic survival model is defined as:
    Surv(x; α, β) = 1-F(x; α, β),
    where α is the LC50 (lethal concentration for 50% of individuals)
    and β is known as the Hill's coefficient, a measure of ultrasensitivity (i.e. how steep is the response curve)
    See [wikipedia on Hill equation](https://en.wikipedia.org/wiki/Hill_equation_(biochemistry))

    Surv(x; α, β) = 1-F(x; α, β) = 1/(1 + (x/ α)^β )
"""
struct SurvLogLogistic{R,W,LC,HC,X} <: SurvRule{R,W}
    "Lethal concentration for 50% of individuals."
    lc50::LC
    "Hill coefficient, or shape of a log logistic function"
    hillcoefficient::HC
    "Exposure: environmental concentration to which the population is exposed to"
    exposure::X
end
function SurvLogLogistic{R,W}(;
    lc50=LC50_PARAM,
    hillcoefficient=HILLCOEFFICIENT_PARAM,
    exposure=EXPOSURE_PARAM,
) where {R,W}
    SurvLogLogistic{R,W}(lc50, hillcoefficient, exposure)
end

@inline function applyrule(data, rule::SurvLogLogistic, population, args...)
    population > zero(population) || return zero(population)
    @fastmath population * ( 1/ (1 + (rule.exposure / rule.lc50)^rule.hillcoefficient) )
end

"""
LogLogistic survival based on an exposure layer
"""
struct SurvLogLogisticMap{R,W,LC,HC,EK,AT} <: SurvMapRule{R,W}
    "Lethal concentration for 50% of individuals."
    lc50::LC
    "Hill coefficient, or shape of a log logistic function"
    hillcoefficient::HC
    "Key for aux data"
    exposurekey::EK
    "Precalculated time interpolation index for aux data"
    auxtimeindex::AT
end
function SurvLogLogisticMap{R,W}(;
    lc50=LC50_PARAM,
    hillcoefficient=HILLCOEFFICIENT_PARAM,
    exposurekey,
    auxtimeindex=1,
) where {R,W}
    SurvLogLogisticMap{R,W}(lc50, hillcoefficient, exposurekey, auxtimeindex)
end

@inline function applyrule(data, rule::SurvLogLogisticMap, population, index, args...)
    population > zero(population) || return zero(population)
    exposure = auxval(data, rule.exposurekey, index..., rule.auxtimeindex) 
    @fastmath population * (1 / (1 + (exposure / rule.lc50)^rule.hillcoefficient))
end

struct GrowthSurvLogLogisticMap3{R,W,LC,HC,CC,GR,EK,AT,TS,S} <: SurvMapRule{R,W}
    "Lethal concentration for 50% of individuals."
    lc50::LC
    "Hill coefficient, or shape of a log logistic function"
    hillcoefficient::HC
    "Carrying capacity for each cell. Not currently scaled by area."
    carrycap::CC
    "Intrinsic rate of growth per timestep"
    intrinsicrate::GR
    "Key for growth rate layer"
    exposurekey::EK
    "Precalculated time interpolation index for aux data"
    auxtimeindex::AT
    "Timestep used in formulation"
    timestep::TS
    "The fractional number of rule timesteps in the current simulation timestep"
    nsteps::S
end
function GrowthSurvLogLogisticMap3{R,W}(;
    lc50=LC50PARAM,
    hillcoefficient=HILLCOEFFICIENT_PARAM,
    carrycap=CARRYCAP_PARAM,
    intrinsicrate=INTRINSICRATE_PARAM,
    exposurekey,
    auxtimeindex=1,
    timestep=nothing,
    nsteps=1.0,
) where {R,W}
    GrowthSurvLogLogisticMap3{R,W}(
        lc50, hillcoefficient, carrycap, intrinsicrate, 
        exposurekey, auxtimeindex, timestep, nsteps,
    )
end

@inline function applyrule(data, rule::GrowthSurvLogLogisticMap3, pops, index, args...)
    growth_carrycap_effect = rule.intrinsicrate * (1 - (sum(pops)) / rule.carrycap)
    exposure = auxval(data, rule.exposurekey, index..., rule.auxtimeindex) 
    map(values(pops), rule.lc50) do pop, lc50
        @fastmath pop * (1 / (1 + (exposure / lc50)^rule.hillcoefficient)) * growth_carrycap_effect
    end
end

# TODO: function named GrowthSurvLogLogisticMap3 is a combination of 2 new function: SurvLogLogisticMapTuple and GrowthMapTuple

struct SurvLogLogisticMapTuple{R,W,LC,HC,EK,AT,TS,S} <: SurvMapRule{R,W}
    "Lethal concentration for 50% of individuals."
    lc50::LC
    "Hill coefficient, or shape of a log logistic function"
    hillcoefficient::HC
    "Key for growth rate layer"
    exposurekey::EK
    "Precalculated time interpolation index for aux data"
    auxtimeindex::AT
    "Timestep used in formulation"
    timestep::TS
    "The fractional number of rule timesteps in the current simulation timestep"
    nsteps::S
end
function SurvLogLogisticMapTuple{R,W}(;
    lc50=LC50PARAM,
    hillcoefficient=HILLCOEFFICIENT_PARAM,
    exposurekey,
    auxtimeindex=1,
    timestep=nothing,
    nsteps=1.0,
) where {R,W}
SurvLogLogisticMapTuple{R,W}(
        lc50, hillcoefficient, exposurekey, auxtimeindex, timestep, nsteps,
    )
end

@inline function applyrule(data, rule::SurvLogLogisticMapTuple, pops, index, args...)
    exposure = auxval(data, rule.exposurekey, index..., rule.auxtimeindex) 
    map(values(pops), rule.lc50) do pop, lc50
        @fastmath pop * (1 / (1 + (exposure / lc50)^rule.hillcoefficient))
    end
end

# ------------------------
struct GrowthMapTuple{R,W,CC,GR,AT,TS,S} <: GrowthMapRule{R,W}
    "Carrying capacity for each cell. Not currently scaled by area."
    carrycap::CC
    "Intrinsic rate of growth per timestep"
    intrinsicrate::GR
    "Precalculated time interpolation index for aux data"
    auxtimeindex::AT
    "Timestep used in formulation"
    timestep::TS
    "The fractional number of rule timesteps in the current simulation timestep"
    nsteps::S
end
function GrowthMapTuple{R,W}(;
    carrycap=CARRYCAP_PARAM,
    intrinsicrate=INTRINSICRATE_PARAM,
    auxtimeindex=1,
    timestep=nothing,
    nsteps=1.0,
) where {R,W}
GrowthMapTuple{R,W}(
        carrycap, intrinsicrate, auxtimeindex, timestep, nsteps,
    )
end

@inline function applyrule(data, rule::GrowthMapTuple, pops, index, args...)
    growth_carrycap_effect = rule.intrinsicrate * (1 - (sum(pops)) / rule.carrycap)
    map(values(pops)) do pop
        @fastmath pop * growth_carrycap_effect
    end
end

"""
LogLogistic survival based on an exposure layer
"""
struct SurvLogLogisticMap2{R,W,HC,EK,LK,AT} <: SurvMapRule{R,W}
    "Hill coefficient, or shape of a log logistic function"
    hillcoefficient::HC
    "Key for exposure data"
    exposurekey::EK
    "Key for growthrate data"
    lc50key::LK
    "Precalculated time interpolation index for aux data"
    auxtimeindex::AT
end
function  SurvLogLogisticMap2{R,W}(;
    hillcoefficient=HILLCOEFFICIENT_PARAM,
    exposurekey,
    lc50key,
    auxtimeindex=1,
) where {R,W}
    SurvLogLogisticMap2{R,W}(hillcoefficient, exposurekey, lc50key, auxtimeindex)
end

@inline function applyrule(data, rule::SurvLogLogisticMap2, population, index, args...)
    population > zero(population) || return zero(population)
    exposure = auxval(data, rule.exposurekey, index..., rule.auxtimeindex) 
    lc50 = auxval(data, rule.lc50key, index..., rule.auxtimeindex) 
    @fastmath population * ( 1/ (1 + (exposure / lc50)^rule.hillcoefficient) )
end

# TODO: Should we gather MaskSurvMap with MaskGrowthMap?
"""
Simple layer mask. Values below a certain threshold are replaced with zero.

"""
struct MaskSurvMap{R,W,ST,EK,AT} <: SurvMapRule{R,W}
    "Minimum viability index."
    threshold::ST
    "Key for growthrate data"
    exposurekey::EK
    "Precalculated time interpolation index for aux data"
    auxtimeindex::AT
end
function MaskSurvMap{R,W}(;
    threshold=Param(0.5; bounds=(0.0, 1.0)),
    exposurekey,
    auxtimeindex=1,
) where {R,W}
    MaskSurvMap{R,W}(threshold, exposurekey, auxtimeindex)
end

@inline function applyrule(data, rule::MaskSurvMap, population, index, args...)
    intrinsicrate = auxval(data, rule.exposurekey, index..., rule.auxtimeindex) 
    intrinsicrate >= rule.threshold ? population : zero(population)
end

const DEGRADATION_RATE = Param(0.5; bounds=(0.0, 100.0))
const PULSE_LEVEL = Param(20.0; bounds=(0.0, 100.0))
const POPULATION_THRESHOLD = Param(20.0; bounds=(0.0, 100.0))
const CROP_OCCUPANCY = Param(5.0; bounds=(0.0, 10.0))

"""
Extends ExposureRule for using heterogenous data.
"""
abstract type Exposure{R,W} <: CellRule{R,W} end

"""
Pulse application
"""
struct Pulsed_Exposure{R,W,PL,TS,S} <: Exposure{R,W}
    "Key for aux layer; pulse of pesticide"
    pulseLevel::PL
    "Timestep used in formulation"
    timestep::TS
    "The fractional number of rule timesteps in the current simulation timestep"
    nsteps::S
end
Pulsed_Exposure{R,W}(;
    pulseLevel=PULSE_LEVEL,
    timestep=nothing,
    nsteps=1.0,
) where {R,W} = Pulsed_Exposure{R,W}(pulseLevel, timestep, nsteps)

precalcrule(rule::Pulsed_Exposure, data) = precalc_timestep(rule, data)

@inline function applyrule(data, rule::Pulsed_Exposure, pesticide, I)
    pulseLevel = get(data, rule.pulseLevel, I...)
    @fastmath pesticide + pulseLevel
end

"""
Degradation of pesticide exposure
"""
struct Degradation_Exposure{R,W,DR,TS,S} <: Exposure{R,W}
    "Degradation Rate for each cell."
    degradationRate::DR
    "Timestep used in formulation"
    timestep::TS
    "The fractional number of rule timesteps in the current simulation timestep"
    nsteps::S
end
Degradation_Exposure{R,W}(;
    degradationRate=DEGRADATION_RATE,
    timestep=nothing,
    nsteps=1.0,
) where {R,W} = Degradation_Exposure{R,W}(degradationRate, timestep, nsteps)

precalcrule(rule::Degradation_Exposure, data) = precalc_timestep(rule, data)

@inline function applyrule(data, rule::Degradation_Exposure, pesticide, I)
    @fastmath pesticide * exp(- rule.degradationRate * rule.nsteps)
end

"""
Threshold application
"""
struct Threshold_Exposure{R,W,CP,PL,PT,DR,TS,S} <: Exposure{R,W}
    "key layer for occurence of plot"
    crop::CP
    "parameter pulse level"
    pulseLevel::PL
    "parameter threshold"
    popThreshold::PT
    "Degradation Rate for each cell."
    degradationRate::DR
    "Timestep used in formulation"
    timestep::TS
    "The fractional number of rule timesteps in the current simulation timestep"
    nsteps::S
end
Threshold_Exposure{R,W}(;
    crop=CROP_OCCUPANCY,
    pulseLevel=PULSE_LEVEL,
    popThreshold=POPULATION_THRESHOLD,
    degradationRate=DEGRADATION_RATE,
    timestep=nothing,
    nsteps=1.0,
) where {R,W} = Threshold_Exposure{R,W}(crop, pulseLevel, popThreshold,degradationRate, timestep, nsteps)

precalcrule(rule::Threshold_Exposure, data) = precalc_timestep(rule, data)

@inline function applyrule(data, rule::Threshold_Exposure, (pesticide, population), I)
    crop = get(data, rule.crop, I...) 
    crop > zero(crop) || return zero(pesticide) 
    # return pesticide only
    if population < rule.popThreshold
        @fastmath pesticide * exp(- rule.degradationRate * rule.nsteps) 
    else 
        @fastmath pesticide * exp(- rule.degradationRate * rule.nsteps) + rule.pulseLevel
    end
end
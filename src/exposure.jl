const DEGRADATION_RATE = Param(0.5; bounds=(0.0, 100.0))
const PULSE_LEVEL = Param(20.0; bounds=(0.0, 100.0))
const POPULATION_THRESHOLD = Param(20.0; bounds=(0.0, 100.0))
const CROP_OCCUPANCY = Param(5.0; bounds=(0.0, 10.0))
const ROTATION_NUMBER = Param(1.0; bounds=(0.0, 1000.0))
const ROTATION_SIZE = Param(1.0; bounds=(0.0, 1000.0))

"""
Extends CellRule for rules of exposure dynamics

For best performance these should be chained with other
CellRule or following an NeighborhoodRule.
"""

abstract type Exposure{R,W} <: CellRule{R,W} end

"""
Exposure by treatment depending on a threshold layer (population size) set with `popThreshold`.
Treatment is applied only on `crop` fields.
The amount of pesticide `pulseLevel` is applied once this threshold is passed.
At every time step the pesticide is degradated at rate `degradationRate`

Read: Tuple(pesticide, population)
Written: Tuple(pesticide)
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
) where {R,W} = Threshold_Exposure{R,W}(crop, pulseLevel, popThreshold, degradationRate, timestep, nsteps)

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

"""
Exposure by several treatments in rotation depending on a threshold layer (population size) set with `popThreshold`.
Treatment is applied only on `crop` fields.
`rotationSize` is the total number of treatments, and `rotationNumber` is the number of treatments.

Example: with 3 treatments, `rotationSize=3` for each treatment rule and `rotationNumber` is set to 1, 2 and 3 respectivelly.

The amount of pesticide `pulseLevel` is applied once the threshold `popThreshold` is passed given by second layer.
At every time step the pesticide is degradated at rate `degradationRate`

Read: Tuple(pesticide, population, rotation)
Written: Tuple(pesticide, rotation)
"""
struct Rotation_Exposure{R,W,CP,RS,RN, PL,PT,DR,TS,S} <: Exposure{R,W}
    "key layer for occurence of plot"
    crop::CP
    "Total number of treatments"
    rotationSize::RS
    "parameter rotation index"
    rotationNumber::RN
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
Rotation_Exposure{R,W}(;
    crop=CROP_OCCUPANCY,
    rotationSize=ROTATION_SIZE,
    rotationNumber=ROTATION_NUMBER,
    pulseLevel=PULSE_LEVEL,
    popThreshold=POPULATION_THRESHOLD,
    degradationRate=DEGRADATION_RATE,
    timestep=nothing,
    nsteps=1.0,
) where {R,W} = Rotation_Exposure{R,W}(crop, rotationSize, rotationNumber, pulseLevel, popThreshold, degradationRate, timestep, nsteps)

precalcrule(rule::Rotation_Exposure, data) = precalc_timestep(rule, data)

@inline function applyrule(data, rule::Rotation_Exposure, (pesticide, population, rotation), I)
    crop = get(data, rule.crop, I...) 
    crop > zero(crop) || return (zero(pesticide), rotation)
    # return pesticide and rotation
    # println("rule.rotationNumber $(rule.rotationNumber) == rotation $(rotation.rotationStep)")
    # println("rotation.rotationStep $(rotation.timeStep) == currentframe $(currentframe(data))")
    # println(
    #     $(population >= rule.popThreshold),
    #     $(mod(rotation.rotationStep, rule.rotationSize) == mod(rule.rotationNumber, rule.rotationSize)),
    #     $(rotation.timeStep != currentframe(data))")
    if  ( population >= rule.popThreshold &&
        mod(rotation.rotationStep, rule.rotationSize) == mod(rule.rotationNumber, rule.rotationSize) &&
        rotation.timeStep != currentframe(data) )
        # println("IN $I")
        @fastmath (
            pesticide * exp(- rule.degradationRate * rule.nsteps) + rule.pulseLevel,
            updateStep(rotation + 1, currentframe(data))
        )
    else
        # println("OUT $I")
        @fastmath (
            pesticide * exp(- rule.degradationRate * rule.nsteps),
            rotation
        )
    end
end

struct RotationStruct{RS,TS}
    rotationStep::RS
    timeStep::TS
end

# All numerical operation are on the first element
Base.:*(rs::RotationStruct, x::Number) = RotationStruct(x * rs.rotationStep, rs.timeStep) 
Base.:*(x::Number, rs::RotationStruct) = RotationStruct(x * rs.rotationStep, rs.timeStep)
Base.:+(rs::RotationStruct, x::Number) = RotationStruct(x + rs.rotationStep, rs.timeStep) 
Base.:+(x::Number, rs::RotationStruct) = RotationStruct(x + rs.rotationStep, rs.timeStep)
Base.:+(rs1::RotationStruct, rs2::RotationStruct) = RotationStruct(rs1.rotationStep + rs2.rotationStep, rs1.timeStep + rs2.timeStep)
Base.:-(rs::RotationStruct, x::Number) = RotationStruct(rs.rotationStep - x, rs.timeStep) 
Base.:-(x::Number, rs::RotationStruct) = RotationStruct(x - rs.rotationStep, rs.timeStep)
Base.:-(rs1::RotationStruct, rs2::RotationStruct) = RotationStruct(rs1.rotationStep - rs2.rotationStep, rs1.timeStep - rs2.timeStep)

Base.mod(rs::RotationStruct, x::Number) = mod(rs.rotationStep, x)

isEqualStep(rs::RotationStruct, y::Number) = rs.timeStep == y ? true : false
isEqualStep(rs::RotationStruct, y::DateTime) = rs.timeStep == y ? true : false
isNotequalStep(rs::RotationStruct, y::Number) = rs.timeStep == y ? false : true
isNotequalStep(rs::RotationStruct, y::DateTime) = rs.timeStep == y ? false : true
isNotequalStep(rs::RotationStruct, y::T) where {T} = rs.timeStep == y ? false : true
updateStep(rs::RotationStruct, y::Number)  = RotationStruct(rs.rotationStep, y)
updateStep(rs::RotationStruct, y::DateTime) = RotationStruct(rs.rotationStep, y)

Base.zero(::Type{<:RotationStruct{T1,T2}}) where {T1,T2} = RotationStruct(zero(T1), zero(T2))

Rotation_InitGrid(tspan_start, dims) = fill(RotationStruct(1,tspan_start), dims)


const DEGRADATION_RATE = Param(0.5; bounds=(0.0, 100.0))
const PULSE_LEVEL = Param(20.0; bounds=(0.0, 100.0))
const POPULATION_THRESHOLD = Param(20.0; bounds=(0.0, 100.0))
const CROP_OCCUPANCY = Param(5.0; bounds=(0.0, 10.0))
const ROTATION_NUMBER = Param(1.0; bounds=(0.0, 1000.0))
const ROTATION_SIZE = Param(1.0; bounds=(0.0, 1000.0))

"""
Abstract supertype of `CellRule` for rules of exposure dynamics
"""

abstract type Exposure{R,W} <: CellRule{R,W} end

"""
    ThresholdExposure{R,W}(;
        crop,  pulselevel,  popthreshold,  degradationrate,  timestep,
    )
Exposure by treatment depending on a threshold layer (population size) set with
`popthreshold`.
Treatment is applied only on `crop` fields (value of `crop` cell > 0).
The amount of treatment `pulselevel` is applied once this threshold is passed.
At every time step the pesticide is degradated at rate `degradationrate`

# Keyword Arguments

- `crop`: Crop field. May be a `Number`, an [`Aux`](@ref) array or another [`Grid`](@ref).
Treatment is applied when `crop > zero(crop) || return zero(pesticide) `.
- `pulselevel`: Pulse level of treatment. May be a `Number`, an [`Aux`](@ref) array or
another [`Grid`](@ref).
- `popthreshold`: Threshold population level triggering the treatment release. May be a
 `Number`, an [`Aux`](@ref) array or another [`Grid`](@ref) of the same type as the
 `population` layer (Read is a Tuple(pesticide, population)) since used in the comparison:
 `population < rule.popthreshold`.
- `degradationrate`: Rate of degradation of the treatment. May be a `Number`, an
 [`Aux`](@ref) array or another [`Grid`](@ref) of a type compatible with `pulselevel`.
- `timestep`: Time step for the exposure, in a type compatible with the simulation `tspan`.

Pass grid `Symbol`s to `R` or both `R` and `W` type parameters to use to specific grids.
Read: Tuple(pesticide, population)
Written: Tuple(pesticide)
"""
struct ThresholdExposure{R,W,CP,PL,PT,DR,TS,S} <: Exposure{R,W}
    crop::CP
    pulselevel::PL
    popthreshold::PT
    degradationrate::DR
    timestep::TS
    nsteps::S
end
ThresholdExposure{R,W}(;
    crop=CROP_OCCUPANCY,
    pulselevel=PULSE_LEVEL,
    popthreshold=POPULATION_THRESHOLD,
    degradationrate=DEGRADATION_RATE,
    timestep,
) where {R,W} = ThresholdExposure{R,W}(crop, pulselevel, popthreshold, degradationrate, timestep, nothing)

precalcrule(rule::ThresholdExposure, data) = precalc_timestep(rule, data)

@inline function applyrule(data, rule::ThresholdExposure, (pesticide, population), I)
    crop = get(data, rule.crop, I...) 
    crop > zero(crop) || return zero(pesticide) 
    pulse = get(data, rule.pulselevel, I...)
    z = get(data, rule.popthreshold, I...)
    rt = get(data, rule.degradationrate, I...) * rule.nsteps
    # return pesticide only
    if population < z
        @fastmath pesticide * exp(- rt) 
    else 
        @fastmath pesticide * exp(- rt) + pulse
    end
end

"""
    RotationExposure{R,W}(;
        crop, rotationsize, rotationindex, popthreshold, degradationrate, timestep,
    )

Exposure by several treatments in rotation depending on a threshold layer (population size) set with `popthreshold`.
Treatment is applied only on `crop` fields.
`rotationsize` is the total number of treatments, and `rotationindex` is the number (index) of a specific treatment.

Example: with 3 treatments, `rotationsize=3` for each treatment rule and `rotationindex` is set to 1, 2 and 3 respectivelly.

The amount of pesticide `pulselevel` is applied once the threshold `popthreshold` is passed given by second layer.
At every time step the pesticide is degradated at rate `degradationrate`

# Keyword Arguments

- `crop`: Crop field. May be a `Number`, an [`Aux`](@ref) array or another [`Grid`](@ref).
Treatment is applied when `crop > zero(crop) || return zero(pesticide) `.
- `rotationsize`: Total number of type in rotation.  May be a `Number`, an [`Aux`](@ref) array or
another [`Grid`](@ref).
- `rotationindex`: Index of the treatment applied. A type compatible with `rotationsize`.
- `pulselevel`: Pulse level of treatment. May be a `Number`, an [`Aux`](@ref) array or
another [`Grid`](@ref).
- `popthreshold`: Threshold population level triggering the treatment release. May be a
 `Number`, an [`Aux`](@ref) array or another [`Grid`](@ref) of the same type as the
 `population` layer (Read is a Tuple(pesticide, population)) since used in the comparison:
 `population < rule.popthreshold`.
- `degradationrate`: Rate of degradation of the treatment. May be a `Number`, an
 [`Aux`](@ref) array or another [`Grid`](@ref) of a type compatible with `pulselevel`.
- `timestep`: Time step for the exposure, in a type compatible with the simulation `tspan`.

Pass grid `Symbol`s to `R` or both `R` and `W` type parameters to use to specific grids.
Read: Tuple(pesticide, population, rotation)
Written: Tuple(pesticide, rotation)
"""
struct RotationExposure{R,W,CP,RS,RN, PL,PT,DR,TS,S} <: Exposure{R,W}
    crop::CP
    rotationsize::RS
    rotationindex::RN
    pulselevel::PL
    popthreshold::PT
    degradationrate::DR
    timestep::TS
    nsteps::S
end
RotationExposure{R,W}(;
    crop=CROP_OCCUPANCY,
    rotationsize=ROTATION_SIZE,
    rotationindex=ROTATION_NUMBER,
    pulselevel=PULSE_LEVEL,
    popthreshold=POPULATION_THRESHOLD,
    degradationrate=DEGRADATION_RATE,
    timestep,
) where {R,W} = RotationExposure{R,W}(crop, rotationsize, rotationindex, pulselevel, popthreshold, degradationrate, timestep, nothing)

precalcrule(rule::RotationExposure, data) = precalc_timestep(rule, data)

@inline function applyrule(data, rule::RotationExposure, (pesticide, population, rotation), I)
    crop = get(data, rule.crop, I...) 
    crop > zero(crop) || return (zero(pesticide), rotation)
    pulse = get(data, rule.pulselevel, I...)
    z = get(data, rule.popthreshold, I...)
    rt = get(data, rule.degradationrate, I...) * rule.nsteps
    if  ( population >= z &&
        mod(_rotationstep(rotation), rule.rotationsize) == mod(rule.rotationindex, rule.rotationsize) &&
        _timestep(rotation) != currentframe(data) )
        # println("IN $I")
        @fastmath (
            pesticide * exp(- rt) + pulse,
            _updatestep(rotation + 1, currentframe(data))
        )
    else
        # println("OUT $I")
        @fastmath (
            pesticide * exp(- rt),
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

Base.zero(::Type{<:RotationStruct{T1,T2}}) where {T1,T2} = RotationStruct(zero(T1), zero(T2))

_updatestep(rs::RotationStruct, y::Number)  = RotationStruct(rs.rotationStep, y)
_updatestep(rs::RotationStruct, y::DateTime) = RotationStruct(rs.rotationStep, y)
_rotationstep(rs::RotationStruct) = rs.rotationStep
_timestep(rs::RotationStruct) = rs.timeStep

InitRotation(tspan_start, dims) = fill(RotationStruct(1,tspan_start), dims)


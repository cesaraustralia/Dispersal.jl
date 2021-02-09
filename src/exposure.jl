const DEGRADATION_RATE = Param(0.5; bounds=(0.0, 100.0))
const PULSE_LEVEL = Param(20.0; bounds=(0.0, 100.0))
const POPULATION_THRESHOLD = Param(20.0; bounds=(0.0, 100.0))
const ACTIVE = Param(5.0; bounds=(0.0, 10.0))
const ROTATION_NUMBER = Param(1.0; bounds=(0.0, 1000.0))
const ROTATION_SIZE = Param(1.0; bounds=(0.0, 1000.0))

"""
Abstract supertype of `CellRule` for rules of exposure dynamics
"""
abstract type Exposure{R,W} <: CellRule{R,W} end

"""
    ExponentialDecrease(; rate, timestep)
    ExponentialDecrease{R}(; rate, timestep)
    ExponentialDecrease{R,W}(; rate, timestep)

Exponential decrease based on a degradation rate data, using exact solution.

# Keyword Arguments

- `rate`: Degradation rate. May be a `Number`, an [`Aux`](@ref) array or another [`Grid`](@ref).
- `timestep`: Time step for the degradation rate, in a type compatible with the simulation `tspan`.

Pass grid `Symbol`s to `R` or both `R` and `W` type parameters to use to specific grids.
"""
struct ExponentialDecrease{R,W,GR,TS,S} <: CellRule{R,W}
    rate::GR
    timestep::TS
    nsteps::S
end
function ExponentialDecrease{R,W}(; 
    rate=DEGRADATION_RATE, 
    timestep, 
) where {R,W}
ExponentialDecrease{R,W}(rate, timestep, nothing)
end

precalcrule(rule::ExponentialDecrease, data) = precalc_timestep(rule, data)

@inline function applyrule(data, rule::ExponentialDecrease, N, I)
    N > zero(N) || return zero(N)
    rt = get(data, rule.rate, I...) * rule.nsteps
    return @fastmath N * exp(- rt)
end

"""
    ThresholdExposure{R,W}(;
        active,  pulselevel,  popthreshold, timestep,
    )
Exposure depending on a threshold layer (population size) set with
`popthreshold`.
Exposure is applied only on `active` fields (value of `active` cell > 0).
The amount of exposure `pulselevel` is applied once this threshold is passed.

# Keyword Arguments

- `active`: Crop field. Exposure is applied when when `active` is above zero.
- `pulselevel`: Pulse level of exposure.
- `popthreshold`: Threshold population level triggering the exposure release. Must be
   of the same type as the `population` layer (Read is a Tuple(pesticide, population))
   since used in the comparison: `population < rule.popthreshold`.
- `timestep`: Time step for the exposure, in a type compatible with the simulation `tspan`.

Arguments `active`, pulselevel` and `popthreshold` may be 
    a `Number`, an [`Aux`](@ref) array or another [`Grid`](@ref).

Pass grid `Symbol`s to `R` or both `R` and `W` type parameters to use to specific grids.
Read: Tuple(exposure, population)
Written: Tuple(exposure)
"""
struct ThresholdExposure{R,W,CP,PL,PT,TS,S} <: Exposure{R,W}
    active::CP
    pulselevel::PL
    popthreshold::PT
    timestep::TS
    nsteps::S
end
ThresholdExposure{R,W}(;
    active=ACTIVE,
    pulselevel=PULSE_LEVEL,
    popthreshold=POPULATION_THRESHOLD,
    timestep,
) where {R,W} = ThresholdExposure{R,W}(active, pulselevel, popthreshold, timestep, nothing)

precalcrule(rule::ThresholdExposure, data) = precalc_timestep(rule, data)

@inline function applyrule(data, rule::ThresholdExposure, (X, N), I)
    active = get(data, rule.active, I...) 
    active > zero(active) || return zero(X) 
    pulse = get(data, rule.pulselevel, I...)
    z = get(data, rule.popthreshold, I...)
    # return pesticide only
    if N < z
        return X
    else 
        return X + pulse
    end
end

"""
    RotationExposure{R,W}(;
        active, rotationsize, rotationindex, popthreshold, timestep,
    )

Several exposures in rotation depending on a threshold layer (population size) set with `popthreshold`.
Exposure is applied only on `active` fields.
`rotationsize` is the total number of exposures, and `rotationindex` is the number (index) of a specific exposure.

Example: with 3 exposures, `rotationsize=3` for each exposure rule and `rotationindex` is set to 1, 2 and 3 respectivelly.

The amount of pesticide `pulselevel` is applied once the threshold `popthreshold` is passed given by second layer.
At every time step the pesticide is degradated at rate `degradationrate`

# Keyword Arguments

- `active`: Crop field. Exposure is applied when `active` is above zero.
- `rotationsize`: Total number of type in rotation.
- `rotationindex`: Index of the exposure applied. A type compatible with `rotationsize`.
- `pulselevel`: Pulse level of exposure. 
- `popthreshold`: Threshold population level triggering the exposure release. May be of
    the same type as the `population` layer (Read is a Tuple(pesticide, population)) 
    since used in the comparison: `population < rule.popthreshold`.
- `timestep`: Time step for the exposure, in a type compatible with the simulation `tspan`.

Arguments `active`, `rotationsize`, `pulselevel` and `popthreshold` may be 
    a `Number`, an [`Aux`](@ref) array or another [`Grid`](@ref).

Pass grid `Symbol`s to `R` or both `R` and `W` type parameters to use to specific grids.
Read: Tuple(pesticide, population, rotation)
Written: Tuple(pesticide, rotation)
"""
struct RotationExposure{R,W,CP,RS,RN, PL,PT,TS,S} <: Exposure{R,W}
    active::CP
    rotationsize::RS
    rotationindex::RN
    pulselevel::PL
    popthreshold::PT
    timestep::TS
    nsteps::S
end
RotationExposure{R,W}(;
    active=ACTIVE,
    rotationsize=ROTATION_SIZE,
    rotationindex=ROTATION_NUMBER,
    pulselevel=PULSE_LEVEL,
    popthreshold=POPULATION_THRESHOLD,
    timestep,
) where {R,W} = RotationExposure{R,W}(
    active, rotationsize, rotationindex, pulselevel, popthreshold, timestep, nothing
)

precalcrule(rule::RotationExposure, data) = precalc_timestep(rule, data)

@inline function applyrule(data, rule::RotationExposure, (X, N, R), I)
    active = get(data, rule.active, I...) 
    active > zero(active) || return (zero(X), R)
    pulse = get(data, rule.pulselevel, I...)
    z = get(data, rule.popthreshold, I...)
    # Test 3 conditions:
    # 1. population should be over the threshold z
    # 2. `rotationstep` should match the modularity of `rotationindex` compared to the number of rotation = `rotationsize` 
    # 3. `datestep` of rotation should be different to currentframe
    if  (N >= z &&
         mod(_rotationstep(R), rule.rotationsize) == mod(rule.rotationindex, rule.rotationsize) &&
         _datestep(R) != currentframe(data))
        # apply pesticide and update datestep and rotation
        return (X + pulse, @fastmath _updaterotation(R, currentframe(data)))
    else
        return (X, R)
    end
end

struct Rotation{RS,TS}
    rotationstep::RS
    datestep::TS
end

_rotationstep(rs::Rotation) = rs.rotationstep
_datestep(rs::Rotation) = rs.datestep

# All numerical operation are on the first element
Base.:*(rs::Rotation, x::Number) = Rotation(x * rs.rotationstep, rs.datestep) 
Base.:*(x::Number, rs::Rotation) = Rotation(x * rs.rotationstep, rs.datestep)
Base.:/(rs::Rotation, x::Number) = Rotation(rs.rotationstep / x, rs.datestep) 
# Keep datestep of first element
Base.:+(rs1::Rotation, rs2::Rotation) = Rotation(rs1.rotationstep + rs2.rotationstep, rs1.datestep)
# Keep datestep of first element
Base.:-(rs1::Rotation, rs2::Rotation) = Rotation(rs1.rotationstep - rs2.rotationstep, rs1.datestep)

Base.isless(rs1::Rotation, rs2::Rotation) = isless(rs1.rotationstep, rs2.rotationstep)
Base.zero(::Type{<:Rotation{T1,T2}}) where {T1,T2} = Rotation(zero(T1), zero(T2))
Base.oneunit(::Type{<:Rotation{T1,T2}}) where {T1,T2} = Rotation(oneunit(T1), oneunit(T2))

_updaterotation(rs::Rotation, y::Int)  = Rotation(rs.rotationstep+1, y)
_updaterotation(rs::Rotation, y::DateTime) = Rotation(rs.rotationstep+1, y)


"""
    initrotation(tspan_start, dims)

Initialise a `rotation` grid fill of `1` at time `tspan_start`. 

# Keyword Arguments

- `tspan_start`: Starting time to apply `RotationExposure` rules.
- `dims`: an `AbstractVector` defining the size of the grid.

"""
initrotation(tspan_start, dims) = fill(Rotation(1,tspan_start), dims)



const CARRYCAP_PARAM      = Param(100000.0; bounds=(0.0, 10.0))
const INTRINSICRATE_PARAM = Param(0.1,      bounds=(0.0, 10.0))
const THRESHOLD_PARAM     = Param(0.5;      bounds=(0.0, 1.0))

"""
Abstract supertype of `CellRule` for growth dynamics rules.
"""
abstract type GrowthRule{R,W} <: CellRule{R,W} end

"""
    ExponentialGrowth(; rate, timestep)
    ExponentialGrowth{R}(; rate, timestep)
    ExponentialGrowth{R,W}(; rate, timestep)

Exponential growth based on a growth rate data, using exact solution.

# Keyword Arguments

- `rate`: Growth rate. May be a `Number`, an [`Aux`](@ref) array or another [`Grid`](@ref).
- `timestep`: Time step for the growth rate, in a type compatible with the simulation `tspan`.

Pass grid `Symbol`s to `R` or both `R` and `W` type parameters to use to specific grids.
"""
struct ExponentialGrowth{R,W,GR,TS,S} <: GrowthRule{R,W}
    rate::GR
    timestep::TS
    nsteps::S
end
function ExponentialGrowth{R,W}(; 
    rate=INTRINSICRATE_PARAM, 
    timestep, 
) where {R,W}
    ExponentialGrowth{R,W}(rate, timestep, nothing)
end

precalcrule(rule::ExponentialGrowth, data) = precalc_timestep(rule, data)

@inline function applyrule(data, rule::ExponentialGrowth, population, I)
    population > zero(population) || return zero(population)
    intrinsicrate = get(data, rule.rate, I...)
    @fastmath population * exp(intrinsicrate * rule.nsteps)
end

"""
    LogisticGrowth(; rate, carrycap, timestep)
    LogisticGrowth{R}(; rate, carrycap, timestep)
    LogisticGrowth{R,W}(; rate, carrycap, timestep)

Logistic growth based on a growth rate layer, using exact solution.

Saturation only applies with positive growth

# Keyword Arguments

- `rate`: Growth rate. May be a `Number`, an [`Aux`](@ref) array or another [`Grid`](@ref).
- `carrycap`: Carrying capacity. May be a `Number`, an [`Aux`](@ref) array or another [`Grid`](@ref).
- `timestep`: Time step for the growth rate, in a type compatible with the simulation `tspan`.

Pass grid `Symbol`s to `R` or both `R` and `W` type parameters to use to specific grids.
"""
struct LogisticGrowth{R,W,GR,CC,TS,S} <: GrowthRule{R,W}
    rate::GR
    carrycap::CC
    timestep::TS
    nsteps::S
end
function LogisticGrowth{R,W}(;
    rate=INTRINSICRATE_PARAM, 
    carrycap=CARRYCAP_PARAM, 
    timestep, 
) where {R,W}
    LogisticGrowth{R,W}(rate, carrycap, timestep, nothing)
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
    ThresholdGrowth(; rate, threshold)
    ThresholdGrowth{G}(; rate, threshold)
    ThresholdGrowth{R,W}(; rate, threshold)

Simple threshold mask. Values below a certain threshold are replaced with zero.

# Keyword Arguments

- `rate`: Growth rate. May be a `Number`, an [`Aux`](@ref) array or another [`Grid`](@ref).
- `threshold`: Minimum viability threshold below which population falls to zero. 
  May be a `Number`, an [`Aux`](@ref) array or another [`Grid`](@ref).

Pass grid `Symbol`s to `R` or both `R` and `W` type parameters to use to specific grids.
"""
struct ThresholdGrowth{R,W,GR,Th} <: GrowthRule{R,W}
    rate::GR
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

precalc_timestep(rule, data) = precalc_timestep(rule.timestep, rule, data)
precalc_timestep(ruletimestep::DatePeriod, rule, data) =
    @set rule.nsteps = currenttimestep(data) / Millisecond(ruletimestep)
precalc_timestep(ruletimestep, rule, data) =
    @set rule.nsteps = timestep(data) / ruletimestep

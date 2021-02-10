
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

Exponential growth of population size N based on an intrinsic growth rate ``r``, using the exact solution between timesteps ``t`` and ``t-1``:

```
N_t = N_{t-1}e^{r t}
```

# Keyword Arguments

- `rate`: Intrinsic growth rate. May be a `Number`, an [`Aux`](@ref) array or another [`Grid`](@ref).
- `timestep`: Time step for the growth rate calculation, in a type compatible with the simulation `tspan`.

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

@inline function applyrule(data, rule::ExponentialGrowth, N, I)
    N > zero(N) || return zero(N)
    rt = get(data, rule.rate, I...) * rule.nsteps

    return @fastmath N * exp(rt)
end

"""
    LogisticGrowth(; rate, carrycap, timestep)
    LogisticGrowth{R}(; rate, carrycap, timestep)
    LogisticGrowth{R,W}(; rate, carrycap, timestep)

Logistic growth rate of population size N based on an intrinsic growth rate ``r`` and carry capacity ``K``, using the exact solution between timesteps ``t+1`` and ``t``:

```
N_{t+1} = (N_t K) / (N_t + (K - N_t) e^{-rt})
```
Saturation only applies with positive growth.

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

@inline function applyrule(data, rule::LogisticGrowth, N, I)
    N > zero(N) || return zero(N)

    rt = get(data, rule.rate, I...) * rule.nsteps 
    k = get(data, rule.carrycap, I...) 

    if rt > zero(rt)
        return @fastmath (N * k) / (N + (k - N) * exp(-rt))
    else
        return @fastmath N * exp(rt)
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

@inline function applyrule(data, rule::ThresholdGrowth, N, I)
    r = get(data, rule.rate, I...)
    t = get(data, rule.threshold, I...)

    return r >= t ? N : zero(N)
end

precalc_timestep(rule, data) = precalc_timestep(rule.timestep, rule, data)
precalc_timestep(rulestep::DatePeriod, rule, data) =
    @set rule.nsteps = currenttimestep(data) / Millisecond(rulestep)
precalc_timestep(rulestep, rule, data) = @set rule.nsteps = timestep(data) / rulestep

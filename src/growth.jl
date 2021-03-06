
const CARRYCAP_PARAM      = Param(100000.0; bounds=(0.0, 10.0))
const INTRINSICRATE_PARAM = Param(0.1,      bounds=(0.0, 10.0))
const THRESHOLD_PARAM     = Param(0.5;      bounds=(0.0, 1.0))

"""
    GrowthRule <: CellRule

Abstract supertype for growth dynamics rules.
"""
abstract type GrowthRule{R,W} <: CellRule{R,W} end

"""
    ExponentialGrowth <: CellRule

    ExponentialGrowth(; rate, timestep, [nsteps_type])
    ExponentialGrowth{R}(; rate, timestep, [nsteps_type])
    ExponentialGrowth{R,W}(; rate, timestep, [nsteps_type])

Exponential growth of population size N based on an intrinsic growth rate ``r``, using the
exact solution between timesteps ``t`` and ``t-1``:

```math
N_t = N_{t-1}e^{r t}
```

# Keywords

- `rate`: Intrinsic growth rate. May be a `Number`, an `Aux` array or another `Grid`.
- `timestep`: Time step for the growth rate calculation, in a type compatible with the 
    simulation `tspan`.
- `nsteps_type`: Specify the floating point type to use when `nsteps` is generated from the
    timestep, if it is required for type-stability or performance. The default is `Float64`.

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
    nsteps_type=Float64,
) where {R,W}
    ExponentialGrowth{R,W}(rate, timestep, zero(nsteps_type))
end

modifyrule(rule::ExponentialGrowth, data) = precalc_nsteps(rule, data)

@inline function applyrule(data, rule::ExponentialGrowth, N, I)
    N > zero(N) || return zero(N)
    rt = get(data, rule.rate, I...) * rule.nsteps

    return @fastmath N * exp(rt)
end

"""
    LogisticGrowth <: GrowthRule

    LogisticGrowth(; rate, carrycap, timestep, [nsteps_type])
    LogisticGrowth{R}(; rate, carrycap, timestep, [nsteps_type])
    LogisticGrowth{R,W}(; rate, carrycap, timestep, [nsteps_type])

Logistic growth rate of population size ``N`` based on an intrinsic growth rate ``r`` and
carry capacity ``K``, using the exact solution between timesteps ``t+1`` and ``t``:

```math
N_{t+1} = (N_t K) / (N_t + (K - N_t) e^{-rt})
```
Saturation only applies with positive growth.

# Keywords

These may be a `Number`, an `Aux` array or another `Grid`:

- `rate`: Intrinsic growth rate.
- `carrycap`: Carrying capacity.
- `timestep`: Time step for the growth rate, in a type compatible with the simulation `tspan`.
- `nsteps_type`: Specify the floating point type to use when `nsteps` is generated from the
    timestep, if it is required for type-stability or performance. The default is `Float64`.

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
    nsteps_type=Float64,
) where {R,W}
    LogisticGrowth{R,W}(rate, carrycap, timestep, zero(nsteps_type))
end

modifyrule(rule::LogisticGrowth, data) = precalc_nsteps(rule, data)

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
    ThresholdGrowth <: CellRule

    ThresholdGrowth(; rate, threshold)
    ThresholdGrowth{R}(; rate, threshold)
    ThresholdGrowth{R,W}(; rate, threshold)

Simple threshold mask. Values below a certain threshold are replaced with zero.

# Keywords

These may be a `Number`, an `Aux` array or another `Grid`.

- `rate`: Intrinsic growth rate.
- `threshold`: Minimum viability threshold below which population falls to zero.

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

precalc_nsteps(rule, data) = precalc_nsteps(rule.timestep, rule, data)
precalc_nsteps(rulestep::DatePeriod, rule, data) =
    @set rule.nsteps = typeof(rule.nsteps)(currenttimestep(data) / Millisecond(rulestep))
precalc_nsteps(rulestep, rule, data) =
    @set rule.nsteps = typeof(rule.nsteps)(timestep(data) / rulestep)

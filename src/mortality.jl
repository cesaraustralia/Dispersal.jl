const HILLCOEFFICIENT = Param(0.1; bounds=(0.0, 100))
const MEDIAN_LOGLOGISTIC = Param(0.1; bounds=(0.0, 100))
const THRESHOLD_EXPOSURE = Param(0.0; bounds=(0.0, 1e9))

"""
    Mortality <: CellRule

Abstract super type for rules of survival effect
For best performance these should be chained with other
CellRule or following an NeighborhoodRule.
"""
abstract type Mortality{R,W} <: CellRule{R,W} end

"""
    LoglogisticMortality <: Mortality

    LoglogisticMortality(; median, hillcoefficient, timestep, [nsteps_type])
    LoglogisticMortality{R}(; median, hillcoefficient, timestep, [nsteps_type])
    LoglogisticMortality{R,W}(; median, hillcoefficient, timestep, [nsteps_type])

Loglogistic mortality based on median, ``α`` and hill coefficient ``β``.

Cumulative function of loglogistic: 
```math
F(x; α, β) =  x^β/(α^β + x^β)
```
where ``α>0`` is the scale and ``β>0`` is the shape.

# Keywords

- `median`: Median of the loglogistic function
- `hillcoefficient`: Hill's coefficient, a measure of ultrasensitivity (i.e. how steep is the response curve).
    May be a `Number`, an `Aux` array or another `Grid`.
- `timestep`: Time step for the mortality rate, in a type compatible with the simulation `tspan`.

Pass grid `Symbol`s to `R` or both `R` and `W` type parameters to use to specific grids.
`R` is a 2 Grids `NamedTuple` like `Tuple{:population,:exposure}` and `W` return only the first grid `:population`.
"""
struct LoglogisticMortality{R,W,MD,HC,TS,S} <: Mortality{R,W}
    median::MD
    hillcoefficient::HC
    timestep::TS
    nsteps::S
end
function LoglogisticMortality{R,W}(;
    median=MEDIAN_LOGLOGISTIC,
    hillcoefficient=HILLCOEFFICIENT,
    timestep,
    nsteps_type=Float64,
) where {R,W}
    LoglogisticMortality{R,W}(median, hillcoefficient, timestep, zero(nsteps_type))
end

modifyrule(rule::LoglogisticMortality, data) = precalc_nsteps(rule, data)

@inline function applyrule(data, rule::LoglogisticMortality, (N, X), I)
    N > zero(N) || return zero(N)
    β = get(data, rule.hillcoefficient, I...)
    α = get(data, rule.median, I...)
    # NOTE: 1 - x^β / (α^β + x^β) = α^β / (α^β + x^β)
    return N * α^β / (α^β + X^β)
end

"""
    ExponentialMortality <: Mortality

    ExponentialMortality(; rate, threshold, timestep, [nsteps_type])
    ExponentialMortality{R}(; rate, threshold, timestep, [nsteps_type])
    ExponentialMortality{R,W}(; rate, threshold, timestep, [nsteps_type])

Exponential mortality based on exposure threshold and mortality rate parameter, using exact solution.

Exponential mortality based on exposure grid ``X``, an exposure threshold parameter ``z``
    and a mortality rate ``r`` using exact solution between time ``t`` and ```t+1``:

```math
N_{t+1} = N_{t}e^{-r t(X-z)}
```

# Keywords

- `rate`: Mortality rate.
- `threshold`: Exposure threshold under which there is no effect.
- `timestep`: Time step for the growth rate, in a type compatible with the simulation `tspan`.
- `nsteps_type`: Specify the floating point type to use when `nsteps` is generated from the
    timestep, if it is required for type-stability or performance. The default is `Float64`.

`rate` and `threshold` can be a `Number`, an `Aux` array or another Grid`.

Pass grid `Symbol`s to `R` or both `R` and `W` type parameters to use to specific grids.
`R` is a 2 Grids `NamedTuple` like `Tuple{:population,:exposure}` and `W` return only the first grid `:population`.
"""
struct ExponentialMortality{R,W,MR,ZX,TS,S} <: GrowthRule{R,W}
    rate::MR
    threshold::ZX
    timestep::TS
    nsteps::S
end
function ExponentialMortality{R,W}(; 
    rate=INTRINSICRATE_PARAM,
    threshold=THRESHOLD_EXPOSURE,
    timestep, 
    nsteps_type=Float64,
) where {R,W}
    ExponentialMortality{R,W}(rate, threshold, timestep, zero(nsteps_type))
end

modifyrule(rule::ExponentialMortality, data) = precalc_nsteps(rule, data)

@inline function applyrule(data, rule::ExponentialMortality, (N,X), I)
    N > zero(N) || return zero(N)
    rt = get(data, rule.rate, I...) * rule.nsteps
    z = get(data, rule.threshold, I...) 
    return N * exp(- rt * max(0, X-z))
end

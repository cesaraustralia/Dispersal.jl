const HILLCOEFFICIENT = Param(0.1; bounds=(0.0, 100))
const MEDIAN_LOGLOGISTIC = Param(0.1; bounds=(0.0, 100))
const THRESHOLD_EXPOSURE = Param(0.0; bounds=(0.0, 1e9))

"""
Extends CellRule for rules of survival effect
For best performance these should be chained with other
CellRule or following an NeighborhoodRule.
"""
abstract type Mortality{R,W} <: CellRule{R,W} end

"""
    LoglogisticMortality(; median, hillcoefficient, timestep)
    LoglogisticMortality{R}(; median, hillcoefficient, timestep)
    LoglogisticMortality{R,W}(; median, hillcoefficient, timestep)

Loglogistic mortality based on median,  ``α`` and hill coefficient ``β``.

Cumulative function of loglogistic: 
```math
F(x; α, β) =  x^β/(α^β + x^β)
```
where α>0 is the scale and β>0 is the shape.

# Keyword Arguments
- `median`: Median of the loglogistic function
- `hillcoefficient`: Hill's coefficient, a measure of ultrasensitivity (i.e. how steep is the response curve).
    May be a `Number`, an [`Aux`](@ref) array or another [`Grid`](@ref).
- `timestep`: Time step for the mortality rate, in a type compatible with the simulation `tspan`.

Pass grid `Symbol`s to `R` or both `R` and `W` type parameters to use to specific grids.
Read: Tuple(population, exposure)
Written: Tuple(population)
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
) where {R,W}
LoglogisticMortality{R,W}(median, hillcoefficient, timestep, nothing)
end

precalcrule(rule::LoglogisticMortality, data) = precalc_timestep(rule, data)

@inline function applyrule(data, rule::LoglogisticMortality, (N, X), I) # RETURN ONLY POPULATION
    N > zero(N) || return zero(N)
    β = get(data, rule.hillcoefficient, I...)
    α = get(data, rule.median, I...)
    # NOTE: 1 - x^β / (α^β + x^β) = α^β / (α^β + x^β)
    @fastmath N * α^β / (α^β + X^β)
end

"""
    LoglogisticMortalityR(; hillcoefficient, timestep)
    LoglogisticMortalityR{R}(; hillcoefficient, timestep)
    LoglogisticMortalityR{R,W}(; hillcoefficient, timestep)

Loglogistic mortality based on a hill coefficient parameter ``β`` and grid-based median ``α``.
Median is set by another grid reflecting change from resistance mechanisms.

Cumulative function of loglogistic: 
```math
F(x; α, β) =  1/(1 + (x/ α)^{-β} ) 
```
where α>0 is the scale and β>0 is the shape.

# Keyword Arguments
- `hillcoefficient`: Hill's coefficient, a measure of ultrasensitivity (i.e. how steep is the response curve).
    May be a `Number`, an [`Aux`](@ref) array or another [`Grid`](@ref).
- `timestep`: Time step for the mortality rate, in a type compatible with the simulation `tspan`.

Pass grid `Symbol`s to `R` or both `R` and `W` type parameters to use to specific grids.
Read: Tuple(population, exposure, median)
Written: Tuple(population)
"""
struct LoglogisticMortalityR{R,W,HC,TS,S} <: Mortality{R,W}
    hillcoefficient::HC
    timestep::TS
    nsteps::S
end
function LoglogisticMortalityR{R,W}(;
    hillcoefficient=HILLCOEFFICIENT,
    timestep,
) where {R,W}
LoglogisticMortalityR{R,W}(hillcoefficient, timestep, nothing)
end

precalcrule(rule::LoglogisticMortalityR, data) = precalc_timestep(rule, data)

@inline function applyrule(data, rule::LoglogisticMortalityR, (N, X, R), I) # RETURN ONLY POPULATION
    N > zero(N) || return zero(N)
    β = get(data, rule.hillcoefficient, I...)
    @fastmath N * R^β / (R^β + X^β)
end


"""
    ExponentialMortality(; rate, threshold, timestep)
    ExponentialMortality{R}(; rate, threshold, timestep)
    ExponentialMortality{R,W}(; rate, threshold, timestep)

Exponential mortality based on exposure threshold and mortality rate parameter, using exact solution.

Exponential mortality based on exposure grid ``X``, an exposure threshold parameter ``z``
    and a mortality rate ``r`` using exact solution between time ``t`` and ```t+1``:

```math
N_{t+1} = N_{t}e^{-r t*(X-z)}
```

# Keyword Arguments

- `rate`: Mortality rate.
- `threshold`: Exposure threshold under which there is no effect.
- `timestep`: Time step for the growth rate, in a type compatible with the simulation `tspan`.

`rate` and `threshold` can be a `Number`, an [`Aux`](@ref) array or another [`Grid`](@ref).

Pass grid `Symbol`s to `R` or both `R` and `W` type parameters to use to specific grids.
Read: Tuple(population, exposure)
Written: Tuple(population)
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
) where {R,W}
    ExponentialMortality{R,W}(rate, threshold, timestep, nothing)
end

precalcrule(rule::ExponentialMortality, data) = precalc_timestep(rule, data)

@inline function applyrule(data, rule::ExponentialMortality, (N,X), I)
    N > zero(N) || return zero(N)
    rt = get(data, rule.rate, I...) * rule.nsteps
    z = get(data, rule.threshold, I...) 
    return @fastmath N * exp(-rt*max(0,X-z))
end

"""
    ExponentialMortalityR(; threshold, timestep)
    ExponentialMortalityR{R}(; threshold, timestep)
    ExponentialMortalityR{R,W}(; threshold, timestep)

Exponential mortality based on exposure grid ``X``, an exposure threshold parameter ``z``
    and a grid-based mortality rate ``r``.
Rate is given by another grid reflecting change according to resistance mechanisms.

```math
N_{t+1} = N_{t}e^{-r t*(X-z)}
```

# Keyword Arguments

- `threshold`: Exposure threshold under which there is no effect.
- `timestep`: Time step for the growth rate, in a type compatible with the simulation `tspan`.

`rate` and `threshold` can be a `Number`, an [`Aux`](@ref) array or another [`Grid`](@ref).

Pass grid `Symbol`s to `R` or both `R` and `W` type parameters to use to specific grids.
Read: Tuple(population, exposure, rate)
Written: Tuple(population)
"""
struct ExponentialMortalityR{R,W,ZX,TS,S} <: GrowthRule{R,W}
    threshold::ZX
    timestep::TS
    nsteps::S
end
function ExponentialMortalityR{R,W}(; 
    threshold=THRESHOLD_EXPOSURE,
    timestep, 
) where {R,W}
    ExponentialMortalityR{R,W}(threshold, timestep, nothing)
end

precalcrule(rule::ExponentialMortalityR, data) = precalc_timestep(rule, data)

@inline function applyrule(data, rule::ExponentialMortalityR, (N,X,R), I)
    N > zero(N) || return zero(N)
    rt = R * rule.nsteps
    z = get(data, rule.threshold, I...) 
    return @fastmath N * exp(-rt*max(0,X-z))
end

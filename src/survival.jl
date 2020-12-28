const LC50_PARAM = Param(0.0; bounds=(0.0, 1e9))
const HILLCOEFFICIENT_PARAM = Param(0.1; bounds=(0.0, 100))

"""
Extends CellRule for rules of survival effect
For best performance these should be chained with other
CellRule or following an NeighborhoodRule.
"""
abstract type Survival{R,W} <: CellRule{R,W} end

"""
Logistic Survival.

Saturation only applies with positive growth
"""
struct LogisticSurvival{R,W,HC,TS,S} <: Survival{R,W}
    "Hill coefficient, or shape of a log logistic function"
    hillcoefficient::HC
    "Timestep used in the formulation"
    timestep::TS
    "The fractional number of rule timesteps in the current simulation timestep"
    nsteps::S
end
function LogisticSurvival{R,W}(;
    hillcoefficient=HILLCOEFFICIENT_PARAM,
    timestep=nothing,
    nsteps=1.0,
) where {R,W}
    LogisticSurvival{R,W}(hillcoefficient,timestep, nsteps)
end

precalcrule(rule::LogisticSurvival, data) = precalc_timestep(rule, data)

@inline function applyrule(data, rule::LogisticSurvival, (population, popPhenotype, exposure), I) # RETURN ONLY POPULATION
    population > zero(population) || return zero(population)
    phenotype = popPhenotype / population
    @fastmath population * ( 1/ (1 + (exposure / phenotype)^rule.hillcoefficient) )
end
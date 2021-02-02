const HILL_COEFFICIENT = Param(0.0; bounds=(0.0, 100.0))
const DEVIATION_PHENOTYPE = Param(0.0; bounds=(0.0, 100.0))
const DOMINANCE_DEGREE = Param(0.0; bounds=(0.0, 100.0))

"""
Extends Resistance for using heterogenous data.
"""
abstract type Resistance{R,W} <: CellRule{R,W} end

"""
Lande Equation application
"""
struct Lande_Resistance{R,W,HC,DP,DD,TS,S} <: Resistance{R,W} 
    "Hill coefficient, or shape of a log logistic function"
    hillcoefficient::HC
    "Deviation of homozygous to the average phenotype"
    deviationPhenotype::DP
    "AdditDegree of dominance between two alleles"
    dominanceDegree::DD
    "Timestep used in formulation"
    timestep::TS
    "The fractional number of rule timesteps in the current simulation timestep"
    nsteps::S
end

Lande_Resistance{R,W}(;
    hillcoefficient = HILL_COEFFICIENT,
    deviationPhenotype = DEVIATION_PHENOTYPE,
    dominanceDegree = DOMINANCE_DEGREE,
    timestep=nothing,
    nsteps=1.0,
) where {R,W} = Lande_Resistance{R,W}(hillcoefficient, deviationPhenotype, dominanceDegree, timestep, nsteps)

precalcrule(rule::Lande_Resistance, data) = precalc_timestep(rule, data)

@inline function applyrule(data, rule::Lande_Resistance, (popAlleleFreq, popPhenotype, population, exposure), I)
    # RETURN Tuple(popAlleleFreq, popAlleleFreq)
    population > zero(population) || return (zero(popAlleleFreq), zero(popPhenotype))
    # Exposure is a layer
    phenotype = popPhenotype / population
    allele_frequency = popAlleleFreq / population
    @fastmath  allele_frequency = allele_frequency +
            allele_frequency * (1-allele_frequency)*(rule.deviationPhenotype + rule.dominanceDegree*(1-2*allele_frequency))*
            rule.hillcoefficient*(exposure/phenotype)^rule.hillcoefficient / ( phenotype *((exposure/phenotype)^rule.hillcoefficient +1))

    @fastmath phenotype = phenotype + 
            2*allele_frequency*(1-allele_frequency) * (rule.deviationPhenotype + rule.dominanceDegree*(1-2*allele_frequency))^2*
          rule.hillcoefficient*(exposure/phenotype)^rule.hillcoefficient / ( phenotype *((exposure/phenotype)^rule.hillcoefficient +1)) *
          (1-rule.dominanceDegree*rule.hillcoefficient*(exposure/phenotype)^rule.hillcoefficient / ( phenotype *((exposure/phenotype)^rule.hillcoefficient +1)))
    
    return (population * allele_frequency, population * phenotype)
end

const ALLELE_FREQUENCY_PARAM = Param(0.5; bounds=(0.0, 1.0))
const DEVIATION_PHENOTYPE_PARAM = Param(0.5; bounds=(-1.0, 1.0))
const DOMINANCE_DEGREE_PARAM = Param(0.5; bounds=(-1.0, 1.0))

abstract type DeltaAlleleFrequencyRule{R,W} <: CellRule{R,W} end

DynamicGrids.precalcrules(rule::DeltaAlleleFrequencyRule, data) = rule

abstract type DeltaAlleleFrequencyMapRule{R,W} <: DeltaAlleleFrequencyRule{R,W} end

struct DeltaAlleleFrequencySurv{R,W,LC,HC,X,AF,DP,DD} <: DeltaAlleleFrequencyRule{R,W}
    "Lethal concentration for 50% of individuals."
    lc50::LC
    "Hill coefficient, or shape of a log logistic function"
    hillcoefficient::HC
    "Exposure: environmental concentration to which the population is exposed to"
    exposure::X
    "Allele frequency"
    allele_frequency::AF
    "Deviation of homozygous to the average phenotype"
    deviation_phenotype::DP
    "AdditDegree of dominance between two alleles"
    dominance_degree::DD
end
function DeltaAlleleFrequencySurv{R,W}(;
    lc50=LC50_PARAM,
    hillcoefficient=HILLCOEFFICIENT_PARAM,
    exposure=EXPOSURE_PARAM,
    allele_frequency=ALLELE_FREQUENCY_PARAM,
    deviation_phenotype=DEVIATION_PHENOTYPE_PARAM,
    dominance_degree=DOMINANCE_DEGREE_PARAM,
) where {R,W}
    DeltaAlleleFrequencySurv{R,W}(
        lc50, hillcoefficient, exposure, allele_frequency, deviation_phenotype, dominance_degree,
    )
end

@inline function applyrule(data, rule::DeltaAlleleFrequencySurv, allele_frequency, index)
    allele_frequency > zero(allele_frequency) || return zero(allele_frequency)
    @fastmath allele_frequency +
        allele_frequency * (1-allele_frequency)*(rule.deviation_phenotype + rule.dominance_degree*(1-2*allele_frequency))*
        rule.hillcoefficient*(rule.exposure/rule.lc50)^rule.hillcoefficient / ( rule.lc50 *((rule.exposure/rule.lc50)^rule.hillcoefficient +1))
end

struct DeltaAlleleFrequencySurvMap{R,W,LC,HC,DP,DD,EK,AT} <: DeltaAlleleFrequencyMapRule{R,W}
    "Lethal concentration for 50% of individuals."
    lc50::LC
    "Hill coefficient, or shape of a log logistic function"
    hillcoefficient::HC
    "Deviation of homozygous to the average phenotype"
    deviation_phenotype::DP
    "AdditDegree of dominance between two alleles"
    dominance_degree::DD
    "Key for growth rate layer"
    exposurekey::EK
    "Precalculated time interpolation index for aux data"
    auxtimeindex::AT
end
function DeltaAlleleFrequencySurvMap{R,W}(;
    lc50=LC50_PARAM,
    hillcoefficient=HILLCOEFFICIENT_PARAM,
    deviation_phenotype=DEVIATION_PHENOTYPE_PARAM,
    dominance_degree=DOMINANCE_DEGREE_PARAM,
    exposurekey,
    auxtimeindex=1,
) where {R,W}
    DeltaAlleleFrequencySurvMap{R,W}(
        lc50, hillcoefficient, deviation_phenotype, dominance_degree, exposurekey, auxtimeindex,
    )
end

function DynamicGrids.precalcrules(rule::DeltaAlleleFrequencySurvMap, data)
    precalc_auxtimeindex(aux(data, rule.exposurekey), rule, data)
end

@inline function applyrule(data, rule::DeltaAlleleFrequencySurvMap, allele_frequency, index)
    allele_frequency > zero(allele_frequency) || return zero(allele_frequency)
    exposure = auxval(data, rule.exposurekey, index..., rule.auxtimeindex) 
    allele_frequency +
        allele_frequency * (1-allele_frequency)*(rule.deviation_phenotype + rule.dominance_degree*(1-2*allele_frequency))*
        rule.hillcoefficient*(exposure/rule.lc50)^rule.hillcoefficient / ( rule.lc50 *((exposure/rule.lc50)^rule.hillcoefficient +1))
end

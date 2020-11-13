
@defkw struct AlleleFrequency{AF}
    "Allele frequency"
    allele_frequency::AF = Param(0.5; bounds=(0.0, 1.0))
end

@defkw struct DeviationPhenotype{DP}
    "Deviation of homozygous to the average phenotype"
    deviation_phenotype::DP Param(0.5; bounds=(-1.0, 1.0))
end

@defkw struct DominanceDegree{DD}
    "AdditDegree of dominance between two alleles"
    dominance_degree::DD = Param(0.5; bounds=(-1.0, 1.0))
end

isprobvec(p::AbstractVector{<:Real}) =
    all(x -> x ≥ zero(x), p) && isapprox(sum(p), one(eltype(p)))


abstract type DeltaAlleleFrequencyRule{R,W} <: CellRule{R,W} end

DynamicGrids.precalcrules(rule::DeltaAlleleFrequencyRule, data) = rule

abstract type DeltaAlleleFrequencyMapRule{R,W} <: DeltaAlleleFrequencyRule{R,W} end

DynamicGrids.precalcrules(rule::DeltaAlleleFrequencyMapRule, data) = begin
    precalclayer(layer(rule, data), rule, data)
end

# @Exposure @LC50 @HillCoefficient @DeviationPhenotype @DominanceDegree 
struct DeltaAlleleFrequencySurv{R,W} <: DeltaAlleleFrequencyRule{R,W} 
end

@inline applyrule(data, rule::DeltaAlleleFrequencySurv, alleleFrequency, args...) = begin
    alleleFrequency > zero(alleleFrequency) || return zero(alleleFrequency)
    @fastmath alleleFrequency + 
        alleleFrequency * (1-alleleFrequency)*(rule.deviation_phenotype + rule.dominance_degree*(1-2*alleleFrequency))*
        rule.hillcoefficient*(rule.exposure/rule.LC50)^rule.hillcoefficient / ( rule.LC50 *((rule.exposure/rule.LC50)^rule.hillcoefficient +1))
end

# @Layers @LC50 @HillCoefficient @DominanceDegree @DeviationPhenotype 
struct DeltaAlleleFrequencySurvMap{R,W} <: DeltaAlleleFrequencyMapRule{R,W} 
    "AdditDegree of dominance between two alleles"
    dominance_degree::DD = Param(0.5; bounds=(-1.0, 1.0))
    "Deviation of homozygous to the average phenotype"
    deviation_phenotype::DP Param(0.5; bounds=(-1.0, 1.0))
end

@inline applyrule(data, rule::DeltaAlleleFrequencySurvMap, allele_frequency, index, args...) = begin
    allele_frequency > zero(allele_frequency) || return zero(allele_frequency)
    exposure = layer(rule, data, index)
    @fastmath allele_frequency +
        allele_frequency * (1-allele_frequency)*(rule.deviation_phenotype + rule.dominance_degree*(1-2*allele_frequency))*
        rule.hillcoefficient*(exposure/rule.LC50)^rule.hillcoefficient / ( rule.LC50 *((exposure/rule.LC50)^rule.hillcoefficient +1))
end

# @Layers @LC50 @HillCoefficient @DeviationPhenotype @DominanceDegree 
struct DeltaAlleleFrequencySurvMap_noFastmath{R,W} <: DeltaAlleleFrequencyMapRule{R,W} 
    "Deviation of homozygous to the average phenotype"
    deviation_phenotype::DP Param(0.5; bounds=(-1.0, 1.0))
    "AdditDegree of dominance between two alleles"
    dominance_degree::DD = Param(0.5; bounds=(-1.0, 1.0))
end

@inline applyrule(data, rule::DeltaAlleleFrequencySurvMap_noFastmath, alleleFrequency, index, args...) = begin
    alleleFrequency > zero(alleleFrequency) || return zero(alleleFrequency)
    exposure = layer(rule, data, index)
    alleleFrequency +
        allele_frequency * (1-allele_frequency)*(rule.deviation_phenotype + rule.dominance_degree*(1-2*allele_frequency))*
        rule.hillcoefficient*(exposure/rule.LC50)^rule.hillcoefficient / ( rule.LC50 *((exposure/rule.LC50)^rule.hillcoefficient +1))
end

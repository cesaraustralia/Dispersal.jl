
@mix @columns struct AlleleFrequency{AF}
    # Field              | Default  | Flat  | Bounds      | Description
    alleleFrequency::AF  | 0.5      | true  | (0.0, 1.0)  | "Allele frequency"
end

@mix @columns struct DeviationPhenotype{DP}
    # Field              | Default  | Flat  | Bounds      | Description
    deviationPhenotype::DP  | 0.5      | true  | (-1.0, 1.0)  | "Deviation of homozygous to the average phenotype"
end

@mix @columns struct DominanceDegree{DD}
    # Field              | Default  | Flat  | Bounds      | Description
    dominanceDegree::DD  | 0.5      | true  | (-1.0, 1.0)  | "AdditDegree of dominance between two alleles"
end

isprobvec(p::AbstractVector{<:Real}) =
    all(x -> x ≥ zero(x), p) && isapprox(sum(p), one(eltype(p)))

# struct AlleleFrequency{T<:AbstractVector{<:Real}}
#     p::T
#     AlleleFrequency{T}(p) where {T<:AbstractVector{<:Real}} = !isprobvec(p) ? throw(ArgumentError("p = $p is not a probability vector.")) : new{T}(p)
# end
#
# AlleleFrequency(1.0)
# AlleleFrequency([0.1, 0.9])
# AlleleFrequency([0.000001, 0.999999])
# AlleleFrequency([1])
# AlleleFrequency([0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.1])

abstract type DeltaAlleleFrequencyRule{R,W} <: CellRule{R,W} end

DynamicGrids.precalcrules(rule::DeltaAlleleFrequencyRule, data) = rule

abstract type DeltaAlleleFrequencyMapRule{R,W} <: DeltaAlleleFrequencyRule{R,W} end

DynamicGrids.precalcrules(rule::DeltaAlleleFrequencyMapRule, data) = begin
    precalclayer(layer(rule, data), rule, data)
end

@Exposure @LC50 @HillCoefficient @DeviationPhenotype @DominanceDegree struct DeltaAlleleFrequencySurv{R,W} <: DeltaAlleleFrequencyRule{R,W} end

@inline applyrule(data, rule::DeltaAlleleFrequencySurv, alleleFrequency, args...) = begin
    alleleFrequency > zero(alleleFrequency) || return zero(alleleFrequency)
    @fastmath alleleFrequency + 
        alleleFrequency * (1-alleleFrequency)*(rule.deviationPhenotype + rule.dominanceDegree*(1-2*alleleFrequency))*
        rule.hillcoefficient*(rule.exposure/rule.LC50)^rule.hillcoefficient / ( rule.LC50 *((rule.exposure/rule.LC50)^rule.hillcoefficient +1))
end

@Layers @LC50 @HillCoefficient @DominanceDegree @DeviationPhenotype struct DeltaAlleleFrequencySurvMap{R,W} <: DeltaAlleleFrequencyMapRule{R,W} end

@inline applyrule(data, rule::DeltaAlleleFrequencySurvMap, alleleFrequency, index, args...) = begin
    alleleFrequency > zero(alleleFrequency) || return zero(alleleFrequency)
    exposure = layer(rule, data, index)
    @fastmath alleleFrequency +
        alleleFrequency * (1-alleleFrequency)*(rule.deviationPhenotype + rule.dominanceDegree*(1-2*alleleFrequency))*
        rule.hillcoefficient*(exposure/rule.LC50)^rule.hillcoefficient / ( rule.LC50 *((exposure/rule.LC50)^rule.hillcoefficient +1))
end

@Layers @LC50 @HillCoefficient @DeviationPhenotype @DominanceDegree struct DeltaAlleleFrequencySurvMap_noFastmath{R,W} <: DeltaAlleleFrequencyMapRule{R,W} end

@inline applyrule(data, rule::DeltaAlleleFrequencySurvMap_noFastmath, alleleFrequency, index, args...) = begin
    alleleFrequency > zero(alleleFrequency) || return zero(alleleFrequency)
    exposure = layer(rule, data, index)
    alleleFrequency +
        alleleFrequency * (1-alleleFrequency)*(rule.deviationPhenotype + rule.dominanceDegree*(1-2*alleleFrequency))*
        rule.hillcoefficient*(exposure/rule.LC50)^rule.hillcoefficient / ( rule.LC50 *((exposure/rule.LC50)^rule.hillcoefficient +1))
end
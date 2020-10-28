@mix @columns struct AdditiveGeneticVariance{VA}
    # Field                      | Default  | Flat  | Bounds      | Description
    additiveGeneticVariance::VA  | 2.0      | true  | (0.0, 1e9)  | "Additive genetic variance"
end

"""
Extends CellRule for rules of selection gradient effect

For best performance these should be chained with other
CellRule or following an NeighborhoodRule.
"""
abstract type SelectionGradientRule{R,W} <: CellRule{R,W} end

"""
Extends SelectionGradientRule for selction rules using an heterogenous
exposure layer.
"""
abstract type SelectionGradientMapRule{R,W} <: SelectionGradientRule{R,W} end

DynamicGrids.precalcrules(rule::SelectionGradientMapRule, data) = begin
    precalclayer(layer(rule, data), rule, data)
end

DynamicGrids.precalcrules(rule::SelectionGradientRule, data) = rule

@Exposure @HillCoefficient @AdditiveGeneticVariance struct SelectionGradientSurv{R,W} <: SelectionGradientRule{R,W} end

@inline applyrule(data, rule::SelectionGradientSurv, LC50, args...) = begin
    LC50 > zero(LC50) || return zero(LC50)
    @fastmath LC50 + rule.additiveGeneticVariance * rule.hillcoefficient*(rule.exposure/LC50)^rule.hillcoefficient / ( LC50 *((rule.exposure/LC50)^rule.hillcoefficient +1))
end

@Layers @HillCoefficient @AdditiveGeneticVariance struct SelectionGradientSurvMap{R,W} <: SelectionGradientMapRule{R,W} end

@inline applyrule(data, rule::SelectionGradientSurvMap, LC50, index, args...) = begin
    LC50 > zero(LC50) || return zero(LC50)
    exposure = layer(rule, data, index)
    @fastmath LC50 + rule.additiveGeneticVariance * rule.hillcoefficient*(exposure/LC50)^rule.hillcoefficient / ( LC50 *((exposure/LC50)^rule.hillcoefficient +1))
end

# ######### ----- Phenotype handling single locus genotype
struct LandeVariable{AF,PH}
    alleleFrequency::AF
    phenotype::PH
end
# const LV = LandeVariable

Base.:*(lv::LandeVariable, x::Number) = LandeVariable(x * lv.alleleFrequency, x * lv.phenotype) 
Base.:*(x::Number, lv::LandeVariable) = LandeVariable(x * lv.alleleFrequency, x * lv.phenotype)
Base.:+(lv::LandeVariable, x::Number) = LandeVariable(x + lv.alleleFrequency, x + lv.phenotype) 
Base.:+(x::Number, lv::LandeVariable) = LandeVariable(x + lv.alleleFrequency, x + lv.phenotype)
Base.:+(lv1::LandeVariable, lv2::LandeVariable) = LandeVariable(lv1.alleleFrequency + lv2.alleleFrequency, lv1.phenotype + lv2.phenotype)
Base.:-(lv::LandeVariable, x::Number) = LandeVariable(x - lv.alleleFrequency, x - lv.phenotype) 
Base.:-(x::Number, lv::LandeVariable) = LandeVariable(x - lv.alleleFrequency, x - lv.phenotype)
Base.:-(lv1::LandeVariable, lv2::LandeVariable) = LandeVariable(lv1.alleleFrequency - lv2.alleleFrequency, lv1.phenotype - lv2.phenotype)

Base.zero(::Type{<:LandeVariable{T1,T2}}) where {T1,T2} = LandeVariable(zero(T1), zero(T2))

@Exposure @HillCoefficient @DominanceDegree @DeviationPhenotype struct SelectionGradient1locusSurv{R,W} <: SelectionGradientRule{R,W} end

@inline applyrule(data, rule::SelectionGradient1locusSurv, landeVariable, args...) = begin
    @fastmath LandeVariable(
        landeVariable.alleleFrequency +
    landeVariable.alleleFrequency * (1-landeVariable.alleleFrequency)*(rule.deviationPhenotype + rule.dominanceDegree*(1-2*landeVariable.alleleFrequency))*
    rule.hillcoefficient*(rule.exposure/landeVariable.phenotype)^rule.hillcoefficient / ( landeVariable.phenotype *((rule.exposure/landeVariable.phenotype)^rule.hillcoefficient +1)),
    landeVariable.phenotype + 
        2*landeVariable.alleleFrequency*(1-landeVariable.alleleFrequency) * (rule.deviationPhenotype + rule.dominanceDegree*(1-2*landeVariable.alleleFrequency))^2*
      rule.hillcoefficient*(rule.exposure/landeVariable.phenotype)^rule.hillcoefficient / ( landeVariable.phenotype *((rule.exposure/landeVariable.phenotype)^rule.hillcoefficient +1)) *
      (1-rule.dominanceDegree*rule.hillcoefficient*(rule.exposure/landeVariable.phenotype)^rule.hillcoefficient / ( landeVariable.phenotype *((rule.exposure/landeVariable.phenotype)^rule.hillcoefficient +1)))
    )
end

@Layers @HillCoefficient @DominanceDegree @DeviationPhenotype struct SelectionGradient1locusSurvMap{R,W} <: SelectionGradientMapRule{R,W} end

@inline applyrule(data, rule::SelectionGradient1locusSurvMap, landeVariable, index, args...) = begin
    exposure = layer(rule, data, index)
    @fastmath LandeVariable(
        landeVariable.alleleFrequency +
    landeVariable.alleleFrequency * (1-landeVariable.alleleFrequency)*(rule.deviationPhenotype + rule.dominanceDegree*(1-2*landeVariable.alleleFrequency))*
    rule.hillcoefficient*(exposure/landeVariable.phenotype)^rule.hillcoefficient / ( landeVariable.phenotype *((exposure/landeVariable.phenotype)^rule.hillcoefficient +1)),
    landeVariable.phenotype + 
        2*landeVariable.alleleFrequency*(1-landeVariable.alleleFrequency) * (rule.deviationPhenotype + rule.dominanceDegree*(1-2*landeVariable.alleleFrequency))^2*
      rule.hillcoefficient*(exposure/landeVariable.phenotype)^rule.hillcoefficient / ( landeVariable.phenotype *((exposure/landeVariable.phenotype)^rule.hillcoefficient +1)) *
      (1-rule.dominanceDegree*rule.hillcoefficient*(exposure/landeVariable.phenotype)^rule.hillcoefficient / ( landeVariable.phenotype *((exposure/landeVariable.phenotype)^rule.hillcoefficient +1)))
    )
end
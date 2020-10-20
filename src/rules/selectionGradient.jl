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
    @fastmath LC50 + rule.additiveGeneticVariance * rule.hillcoefficient*(rule.exposure/LC50)^rule.hillcoefficient / ( LC50 *((rule.exposure/LC50)^rule.hillcoefficient +1)^2)
end

@Layers @HillCoefficient @AdditiveGeneticVariance struct SelectionGradientSurvMap{R,W} <: SelectionGradientMapRule{R,W} end

@inline applyrule(data, rule::SelectionGradientSurvMap, LC50, index, args...) = begin
    LC50 > zero(LC50) || return zero(LC50)
    exposure = layer(rule, data, index)
    @fastmath LC50 + rule.additiveGeneticVariance * rule.hillcoefficient*(exposure/LC50)^rule.hillcoefficient / ( LC50 *((exposure/LC50)^rule.hillcoefficient +1)^2)
end
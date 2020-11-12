@mix @columns struct LC50{LC}
    # Field           | Default  | Flat  | Bounds      | Description
    LC50::LC          | 10.0     | true  | (0.0, 1e9)  | "Lethal concentration for 50% of individuals."
end

@mix @columns struct HillCoefficient{HC}
    # Field                 | Default  | Flat  | Bounds      | Description
    hillcoefficient::HC     | 0.1      | true  | (0.0, 1e9) | "Hill coefficient, or shape of a log logistic function"
end

@mix @columns struct Exposure{XP}
    # Field           | Default  | Flat  | Bounds      | Description
    exposure::XP      | 100.0      | true  | (0.0, 1e9) | "Exposure: environmental concentration to which the population is exposed to"
end

"""
Extends CellRule for rules of survival effect

For best performance these should be chained with other
CellRule or following an NeighborhoodRule.
"""
abstract type SurvRule{R,W} <: CellRule{R,W} end

"""
Extends SurvRule for survival rules using a heterogenous
concentration layer.
"""

abstract type SurvMapRule{R,W} <: SurvRule{R,W} end

# Base.key(rule::SurvMapRule{R,W,K}) where {R,W,K} = K
# layer(rule::SurvMapRule, data) = axkey(rule)

DynamicGrids.precalcrules(rule::SurvMapRule, data) = begin
    precalclayer(layer(rule, data), rule, data)
end

DynamicGrids.precalcrules(rule::SurvRule, data) = rule

# Exact survival solutions

"""
Simple fixed survival loglogistic using exact solution.

Cumulative function of loglogistic: F(x; α, β) =  1/(1 + (x/ α)^{-β} ),
    where α>0 is the scale and β>0 is the shape.

Dose-response log-logistic survival model is defined as:
    Surv(x; α, β) = 1-F(x; α, β),
    where α is the LC50 (lethal concentration for 50% of individuals)
    and β is known as the Hill's coefficient, a measure of ultrasensitivity (i.e. how steep is the response curve)
    See [wikipedia on Hill equation](https://en.wikipedia.org/wiki/Hill_equation_(biochemistry))    

    Surv(x; α, β) = 1-F(x; α, β) = 1/(1 + (x/ α)^β )

$(FIELDDOCTABLE)
"""
@Exposure @LC50 @HillCoefficient struct SurvLogLogistic{R,W} <: SurvRule{R,W} end

@inline applyrule(data, rule::SurvLogLogistic, population, args...) = begin
    population > zero(population) || return zero(population)
    @fastmath population * ( 1/ (1 + (rule.exposure / rule.LC50)^rule.hillcoefficient) )
end

"""
LogLogistic survival based on an exposure layer
$(FIELDDOCTABLE)
"""
@Layers @LC50 @HillCoefficient struct SurvLogLogisticMap{R,W} <: SurvMapRule{R,W} end

@inline applyrule(data, rule::SurvLogLogisticMap, population, index, args...) = begin
    population > zero(population) || return zero(population)
    exposure = layer(rule, data, index)
    @fastmath population * ( 1/ (1 + (exposure / rule.LC50)^rule.hillcoefficient) )
end

@mix @columns struct LC503{LC1,LC2,LC3}
    # Field           | Default  | Flat  | Bounds      | Description
    LC501::LC1        | 10.0     | true  | (0.0, 1e9)  | "Lethal concentration for 50% of individuals."
    LC502::LC2        | 10.0     | true  | (0.0, 1e9)  | "Lethal concentration for 50% of individuals."
    LC503::LC3        | 10.0     | true  | (0.0, 1e9)  | "Lethal concentration for 50% of individuals."
end

@Layers @LC503 @HillCoefficient @CarryCap @InstrinsicGrowthRate struct GrowthSurvLogLogisticMap3{R,W} <: SurvMapRule{R,W} end

@inline applyrule(data, rule::GrowthSurvLogLogisticMap3, (population1,population2,population3), index, args...) = begin
    growthCarrycapEffect =  rule.intrinsicrate * (1-(population1+population2+population3)/rule.carrycap)
    exposure = layer(rule, data, index)
    @fastmath (
        population1 * ( 1/ (1 + (exposure / rule.LC501)^rule.hillcoefficient) ) * growthCarrycapEffect,
        population2 * ( 1/ (1 + (exposure / rule.LC502)^rule.hillcoefficient) ) * growthCarrycapEffect,
        population3 * ( 1/ (1 + (exposure / rule.LC503)^rule.hillcoefficient) ) * growthCarrycapEffect
    )
end


"""
LogLogistic survival based on an exposure layer
$(FIELDDOCTABLE)
"""
@Layers @HillCoefficient struct SurvLogLogisticMap2{R,W} <: SurvMapRule{R,W} end

@inline applyrule(data, rule::SurvLogLogisticMap2, population, index, args...) = begin
    population > zero(population) || return zero(population)
    exposure = layer(rule, data, index,1)
    LC50 = layer(rule, data, index,2)
    @fastmath population * ( 1/ (1 + (exposure / LC50)^rule.hillcoefficient) )
end


"""
Simple layer mask. Values below a certain threshold are replaced with zero.
$(FIELDDOCTABLE)
"""

#? Should we gather MaskSurvMap with MaskGrowthMap?
@Layers @Timestep struct MaskSurvMap{R,W,ST} <: SurvMapRule{R,W}
    threshold::ST | 0.5 | true  | (0.0, 1.0) | "Minimum viability index."
end

@inline applyrule(data, rule::MaskSurvMap, population, index, args...) =
    layer(rule, data, index) >= rule.threshold ? population : zero(population)
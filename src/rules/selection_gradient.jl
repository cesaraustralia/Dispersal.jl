const ADDITIVE_GV_PARAM = Param(2.0; bounds=(0.0, 1e9))

"""
Extends CellRule for rules of selection gradient effect

For best performance these should be chained with other
CellRule or following an NeighborhoodRule.
"""
abstract type SelectionGradientRule{R,W} <: CellRule{R,W} end

"""
Extends SelectionGradientRule for selction rules using an heterogenous exposure layer.
"""
abstract type SelectionGradientMapRule{R,W} <: SelectionGradientRule{R,W} end

function DynamicGrids.precalcrules(rule::SelectionGradientMapRule, data)
    precalc_auxtimeindex(aux(data, rule.exposurekey), rule, data)
end

DynamicGrids.precalcrules(rule::SelectionGradientRule, data) = rule

struct SelectionGradientSurv{R,W,HC,X,VA} <: SelectionGradientRule{R,W} 
    "Hill coefficient, or shape of a log logistic function"
    hillcoefficient::HC
    "Exposure: environmental concentration to which the population is exposed to"
    exposure::X
    "Additive genetic variance"
    additive_genetic_variance::VA
end
function SelectionGradientSurv{R,W}(;
    hillcoefficient=HILLCOEFFICIENT_PARAM,
    exposure=EXPOSURE_PARAM,
    additive_genetic_variance=ADDITIVE_GV_PARAM,
) where {R,W}
    SelectionGradientSurv{R,W}(hillcoefficient, exposure, additive_genetic_variance)
end

@inline function applyrule(data, rule::SelectionGradientSurv, LC50, index)
    LC50 > zero(LC50) || return zero(LC50)
    @fastmath LC50 + rule.additive_genetic_variance * rule.hillcoefficient*(rule.exposure/LC50)^rule.hillcoefficient / ( LC50 *((rule.exposure/LC50)^rule.hillcoefficient +1))
end

struct SelectionGradientSurvMap{R,W,HC,VA,EK,AT} <: SelectionGradientMapRule{R,W} 
    "Hill coefficient, or shape of a log logistic function"
    hillcoefficient::HC
    "Additive genetic variance"
    additive_genetic_variance::VA
    "Key for exposure layer"
    exposurekey::EK
    "Precalculated time interpolation index for aux data"
    auxtimeindex::AT
end
function SelectionGradientSurvMap{R,W}(;
    hillcoefficient=HILLCOEFFICIENT_PARAM,
    additive_genetic_variance=ADDITIVE_GV_PARAM,
    exposurekey,
    auxtimeindex=1,
) where {R,W}
    SelectionGradientSurvMap{R,W}(
        hillcoefficient, additive_genetic_variance, exposurekey, auxtimeindex
    )
end

@inline function applyrule(data, rule::SelectionGradientSurvMap, LC50, index)
    LC50 > zero(LC50) || return zero(LC50)
    exposure = auxval(data, rule.exposurekey, index..., rule.auxtimeindex) 
    @fastmath LC50 + rule.additive_genetic_variance * rule.hillcoefficient*(exposure/LC50)^rule.hillcoefficient / ( LC50 *((exposure/LC50)^rule.hillcoefficient +1))
end

# ######### ----- Phenotype handling single locus genotype
struct LandeVariable{AF,PH}
    allele_frequency::AF
    phenotype::PH
end
# const LV = LandeVariable

Base.:*(lv::LandeVariable, x::Number) = LandeVariable(x * lv.allele_frequency, x * lv.phenotype) 
Base.:*(x::Number, lv::LandeVariable) = LandeVariable(x * lv.allele_frequency, x * lv.phenotype)
Base.:+(lv::LandeVariable, x::Number) = LandeVariable(x + lv.allele_frequency, x + lv.phenotype) 
Base.:+(x::Number, lv::LandeVariable) = LandeVariable(x + lv.allele_frequency, x + lv.phenotype)
Base.:+(lv1::LandeVariable, lv2::LandeVariable) = LandeVariable(lv1.allele_frequency + lv2.allele_frequency, lv1.phenotype + lv2.phenotype)
Base.:-(lv::LandeVariable, x::Number) = LandeVariable(x - lv.allele_frequency, x - lv.phenotype) 
Base.:-(x::Number, lv::LandeVariable) = LandeVariable(x - lv.allele_frequency, x - lv.phenotype)
Base.:-(lv1::LandeVariable, lv2::LandeVariable) = LandeVariable(lv1.allele_frequency - lv2.allele_frequency, lv1.phenotype - lv2.phenotype)

Base.zero(::Type{<:LandeVariable{T1,T2}}) where {T1,T2} = LandeVariable(zero(T1), zero(T2))

struct SelectionGradient1locusSurv{R,W,HC,X,DP,DD} <: SelectionGradientRule{R,W} 
    "Hill coefficient, or shape of a log logistic function"
    hillcoefficient::HC
    "Exposure: environmental concentration to which the population is exposed to"
    exposure::X
    "Deviation of homozygous to the average phenotype"
    deviation_phenotype::DP
    "AdditDegree of dominance between two alleles"
    dominance_degree::DD
end
function SelectionGradient1locusSurv{R,W}(;
    hillcoefficient=HILLCOEFFICIENT_PARAM,
    exposure=EXPOSURE_PARAM,
    dominance_degree=DOMINANCE_DEGREE_PARAM,
    deviation_phenotype=DEVIATION_PHENOTYPE_PARAM,
) where {R,W}
    SelectionGradient1locusSurv{R,W}(hillcoefficient, exposure, deviation_phenotype, dominance_degree)
end

@inline function applyrule(data, rule::SelectionGradient1locusSurv, lande_variable, index)
    @fastmath LandeVariable(
        lande_variable.allele_frequency +
        lande_variable.allele_frequency * (1-lande_variable.allele_frequency)*(rule.deviation_phenotype + rule.dominance_degree*(1-2*lande_variable.allele_frequency))*
        rule.hillcoefficient*(rule.exposure/lande_variable.phenotype)^rule.hillcoefficient / ( lande_variable.phenotype *((rule.exposure/lande_variable.phenotype)^rule.hillcoefficient +1))
        ,
        lande_variable.phenotype + 
            2*lande_variable.allele_frequency*(1-lande_variable.allele_frequency) * (rule.deviation_phenotype + rule.dominance_degree*(1-2*lande_variable.allele_frequency))^2*
          rule.hillcoefficient*(rule.exposure/lande_variable.phenotype)^rule.hillcoefficient / ( lande_variable.phenotype *((rule.exposure/lande_variable.phenotype)^rule.hillcoefficient +1)) *
          (1-rule.dominance_degree*rule.hillcoefficient*(rule.exposure/lande_variable.phenotype)^rule.hillcoefficient / ( lande_variable.phenotype *((rule.exposure/lande_variable.phenotype)^rule.hillcoefficient +1)))
    )
end

struct SelectionGradient1locusSurvMap{R,W,HC,DP,DD,EK,AT} <: SelectionGradientMapRule{R,W} 
    "Hill coefficient, or shape of a log logistic function"
    hillcoefficient::HC
    "Deviation of homozygous to the average phenotype"
    deviation_phenotype::DP
    "AdditDegree of dominance between two alleles"
    dominance_degree::DD
    "Key for aux data"
    exposurekey::EK
    "Precalculated time interpolation index for aux data"
    auxtimeindex::AT
end
function SelectionGradient1locusSurvMap{R,W}(;
    hillcoefficient=HILLCOEFFICIENT_PARAM,
    deviation_phenotype=DEVIATION_PHENOTYPE_PARAM,
    dominance_degree=DOMINANCE_DEGREE_PARAM,
    exposurekey,
    auxtimeindex=1,
) where {R,W}
    SelectionGradient1locusSurvMap{R,W}(
        hillcoefficient, deviation_phenotype, dominance_degree, exposurekey, auxtimeindex
    )
end

@inline function applyrule(data, rule::SelectionGradient1locusSurvMap, lande_variable, index)
    exposure = auxval(data, rule.exposurekey, index..., rule.auxtimeindex) 
    @fastmath LandeVariable(
        lande_variable.allele_frequency +
        lande_variable.allele_frequency * (1-lande_variable.allele_frequency)*(rule.deviation_phenotype + rule.dominance_degree*(1-2*lande_variable.allele_frequency))*
        rule.hillcoefficient*(exposure/lande_variable.phenotype)^rule.hillcoefficient / ( lande_variable.phenotype *((exposure/lande_variable.phenotype)^rule.hillcoefficient +1)),

        lande_variable.phenotype + 
            2*lande_variable.allele_frequency*(1-lande_variable.allele_frequency) * (rule.deviation_phenotype + rule.dominance_degree*(1-2*lande_variable.allele_frequency))^2*
          rule.hillcoefficient*(exposure/lande_variable.phenotype)^rule.hillcoefficient / ( lande_variable.phenotype *((exposure/lande_variable.phenotype)^rule.hillcoefficient +1)) *
          (1-rule.dominance_degree*rule.hillcoefficient*(exposure/lande_variable.phenotype)^rule.hillcoefficient / ( lande_variable.phenotype *((exposure/lande_variable.phenotype)^rule.hillcoefficient +1)))
    )
end



struct SelectionGradientMapTuple{R,W,HC,DP,DD,EK,AT} <: SelectionGradientMapRule{R,W} 
    "Hill coefficient, or shape of a log logistic function"
    hillcoefficient::HC
    "Deviation of homozygous to the average phenotype"
    deviation_phenotype::DP
    "AdditDegree of dominance between two alleles"
    dominance_degree::DD
    "Key for aux data"
    exposurekey::EK
    "Precalculated time interpolation index for aux data"
    auxtimeindex::AT
end
function SelectionGradientMapTuple{R,W}(;
    hillcoefficient=HILLCOEFFICIENT_PARAM,
    deviation_phenotype=DEVIATION_PHENOTYPE_PARAM,
    dominance_degree=DOMINANCE_DEGREE_PARAM,
    exposurekey,
    auxtimeindex=1,
) where {R,W}
SelectionGradientMapTuple{R,W}(
        hillcoefficient, deviation_phenotype, dominance_degree, exposurekey, auxtimeindex
    )
end

@inline function applyrule(data, rule::SelectionGradientMapTuple, (allele_frequency, phenotype), index)
    exposure = auxval(data, rule.exposurekey, index..., rule.auxtimeindex) 
    @fastmath (allele_frequency +
        allele_frequency * (1-allele_frequency)*(rule.deviation_phenotype + rule.dominance_degree*(1-2*allele_frequency))*
        rule.hillcoefficient*(exposure/phenotype)^rule.hillcoefficient / ( phenotype *((exposure/phenotype)^rule.hillcoefficient +1)),

        phenotype + 
            2*allele_frequency*(1-allele_frequency) * (rule.deviation_phenotype + rule.dominance_degree*(1-2*allele_frequency))^2*
          rule.hillcoefficient*(exposure/phenotype)^rule.hillcoefficient / ( phenotype *((exposure/phenotype)^rule.hillcoefficient +1)) *
          (1-rule.dominance_degree*rule.hillcoefficient*(exposure/phenotype)^rule.hillcoefficient / ( phenotype *((exposure/phenotype)^rule.hillcoefficient +1)))
    )
end

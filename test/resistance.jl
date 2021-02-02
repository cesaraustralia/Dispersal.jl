using Dispersal
using DynamicGrids, Unitful, Dates
using DimensionalData
using Test

@testset "Test resistance Lande equations" begin

population = [  100000.0 100.0 0.0;
                0.0 1000.0 0.0;
                0.0 0.0 100000.0]

pFreqinit = [1.0 0.0 0.0;
             0.0 0.0 0.0;
             0.0 0.0 0.0]
PopAlleleFreq = population .* pFreqinit # all are SS

meanPhenotype = 30.0
deviationPhenotype = 5.0
dominanceDegree = 10.0
PopPhenotype = population .* (meanPhenotype .+ pFreqinit.^2 .* deviationPhenotype + 2 .* pFreqinit .*(1 .-pFreqinit) .* dominanceDegree - (1 .-pFreqinit).^2 .* deviationPhenotype)

init = (
    pesticide =  [0.0 0.0 0.0; 
                0.0 0.0 0.0;
                0.0 0.0 0.0],
    population = population,
    popAlleleFreq = PopAlleleFreq,
    popPhenotype = PopPhenotype,
                )

crop = [0.0 1.0 1.0; 
        0.0 10.0 0.0;
        0.0 0.0 10.0]
    
hillcoefficient = 0.5
output = ArrayOutput(init; tspan=1:3,  aux=(crop=crop,))
ruleTreatment = Threshold_Exposure{Tuple{:pesticide,:population}, :pesticide}(crop=Aux(:crop), pulseLevel=20, popThreshold = 1000)
ruleGrowth = LogisticGrowth{:population,:population}(rate=5)
ruleSurvival = LogisticSurvival{Tuple{:population,:popPhenotype,:pesticide}, :population}(hillcoefficient=hillcoefficient)
ruleResistance = Lande_Resistance{Tuple{:popAlleleFreq,:popPhenotype,:population,:pesticide}, Tuple{:popAlleleFreq,:popPhenotype}}(
    hillcoefficient=hillcoefficient, deviationPhenotype=deviationPhenotype, dominanceDegree=dominanceDegree)

rule = Ruleset(ruleGrowth, ruleTreatment, ruleSurvival, ruleResistance)
sim!(output, rule)
output[Ti(1)]
output[Ti(2)]
output[Ti(3)]

end


population = [100000.0 100.0 0.0;
              0.0 1000.0 0.0;
              0.0 0.0 100000.0]

pFreqinit = [1.0 0.0 0.0;
             0.0 0.0 0.0;
             0.0 0.0 0.0]
PopAlleleFreq = population .* pFreqinit # all are SS

meanPhenotype = 30.0
deviationPhenotype = 5.0
dominanceDegree = 10.0
PopPhenotype = population .* (meanPhenotype .+ pFreqinit.^2 .* deviationPhenotype + 2 .* pFreqinit .*(1 .-pFreqinit) .* dominanceDegree - (1 .-pFreqinit).^2 .* deviationPhenotype)

init = (
    pesticide =  [0.0 0.0 0.0; 
                0.0 0.0 0.0;
                0.0 0.0 0.0],
    population = population,
    popAlleleFreq = PopAlleleFreq,
    popPhenotype = PopPhenotype,
                )

crop = [0.0 1.0 1.0; 
        0.0 10.0 0.0;
        0.0 0.0 10.0]
    
hillcoefficient = 0.5
output = ArrayOutput(init; tspan=1:3,  aux=(crop=crop,))
ruleTreatment = Threshold_Exposure{Tuple{:pesticide,:population}, :pesticide}(crop=Aux(:crop), pulseLevel=20, popThreshold = 1000)
ruleSurvival = LogisticSurvival{Tuple{:population,:popPhenotype,:pesticide}, :population}(hillcoefficient=hillcoefficient)

rule = Ruleset(ruleTreatment, ruleSurvival)
sim!(output, rule)
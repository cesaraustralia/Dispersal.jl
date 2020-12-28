using Dispersal
using DynamicGrids, Unitful, Dates
using DimensionalData
using Test

@testset "Test generate exposure stack" begin

    pesticide =  [0.0 0.0 0.0; 
                  0.0 0.0 0.0;
                  0.0 0.0 0.0]

    pulse = [10.0 10.0 0.0;
                0.0 1.0 0.0;
                0.0 0.0 10.0]

    output = ArrayOutput(pesticide; tspan=1:1:6, aux=(pulse=pulse,))
    
    rulePulse = Pulsed_Exposure(pulseLevel=Aux(:pulse))
    ruleDegradation = Degradation_Exposure(degradationRate=0.2)

    rule = Ruleset(Chain(rulePulse,ruleDegradation))
    sim!(output, rule)
    output[Ti(1)]
    output[Ti(2)]
    output[Ti(3)]
    output[Ti(4)]

    output = ArrayOutput(pesticide; tspan=1:1:6, aux=(pulse=pulse,))
    rule =  Ruleset(Pulsed_Exposure(pulseLevel=Aux(:pulse), timestep = 1:2:6))
    sim!(output, rule)
    output[Ti(1)]
    output[Ti(2)]
    output[Ti(3)]
    output[Ti(4)]
    output[Ti(5)]
 
end    


@testset "Test exposure population" begin

    init = (
        pesticide =  [0.0 0.0 0.0; 
                    0.0 0.0 0.0;
                    0.0 0.0 0.0],
        population = [100000.0 100.0 0.0;
                    0.0 1000.0 0.0;
                    0.0 0.0 100000.0], )

    crop = [0.0 1.0 1.0; 
    0.0 10.0 0.0;
    0.0 0.0 10.0]
        
    output = ArrayOutput(init; tspan=1:3,  aux=(crop=crop,))
    ruleTreatment = Threshold_Exposure{Tuple{:pesticide,:population}, :pesticide}(crop=Aux(:crop), pulseLevel=20, popThreshold = 1000)
    rule = Ruleset(ruleTreatment)
    sim!(output, rule)
    output[Ti(1)]
    output[Ti(2)]
    output[Ti(3)]

    # If cropping change allong time
    crop = reshape( # TODO: something that change with time
        repeat([0.0 1.0 1.0; 
                0.0 10.0 0.0;
                0.0 0.0 10.0], outer=(1,3) ), (3,3,3))
    cropDD = DimArray(crop, (X, Y, Ti));

    output = ArrayOutput(init; tspan=1:3,  aux=(cropDD=cropDD,))
    
    ruleTreatment = Threshold_Exposure{Tuple{:pesticide,:population}, :pesticide}(crop=Aux(:cropDD), pulseLevel=20, popThreshold = 1000)

    rule = Ruleset(ruleTreatment)
    sim!(output, rule)
    output[Ti(1)].pesticide ; output[Ti(1)].population
    output[Ti(2)].pesticide
    output[Ti(3)].pesticide

    
    output = ArrayOutput(init; tspan=1:3,  aux=(cropDD=cropDD,))
    
    ruleTreatment = Threshold_Exposure{Tuple{:pesticide,:population}, :pesticide}(crop=Aux(:cropDD), pulseLevel=20, popThreshold = 1000)
    ruleGrowth = LogisticGrowth{:population,:population}(rate=5)
    rule = Ruleset(ruleGrowth,ruleTreatment)
    sim!(output, rule)
    output[Ti(1)]
    output[Ti(2)]
    output[Ti(3)]

end   
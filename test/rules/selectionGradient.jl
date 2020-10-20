using Dispersal, Test

@testset "selection gradient on LC50" begin
    init =  [10.0 40.0 70.0;
             20.0 50.0 80.0;
             30.0 60.0 90.0]

    output = ArrayOutput(init; tspan=1:3)
    rule = Ruleset(SelectionGradientSurv(exposure = 10, hillcoefficient=2.0, additiveGeneticVariance = 1.0))
    sim!(output, rule)

    @test output[1] == test1
    @test output[2] ≈ test2 atol=1e-4
    @test output[3] ≈ test3 atol=1e-4

end

@testset "selection gradient on LC50 Map" begin
    LC50init =  [10.0 40.0 70.0;
                20.0 50.0 80.0;
                30.0 60.0 90.0]

    exposure = [50.0 50.0 50.0;
                20.0 20.0 20.0;
                0.0 0.0 0.0]

    output = ArrayOutput(LC50init; tspan=1:3, aux=(exposure=exposure,))
    rule = Ruleset(SelectionGradientSurvMap(layerkey=Val(:exposure), hillcoefficient=2.0, additiveGeneticVariance = 1.0))
    sim!(output, rule)

    @test output[1] == test1
    @test output[2] ≈ [10.0074  40.0119  70.0064;
                       20.025   50.0048  80.0014;
                       30.0     60.0     90.0] atol=1e-4
    @test output[3] ≈ [10.0148  40.0238  70.0128;
                       20.05    50.0095  80.0028;
                       30.0     60.0     90.0] atol=1e-4
end
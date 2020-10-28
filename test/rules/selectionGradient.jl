using Dispersal, Test

@testset "selection gradient on LC50" begin
    init =  [10.0 40.0 70.0;
             20.0 50.0 80.0;
             30.0 60.0 90.0]

    output = ArrayOutput(init; tspan=1:3)
    rule = Ruleset(SelectionGradientSurv(exposure = 10, hillcoefficient=2.0, additiveGeneticVariance = 1.0))
    sim!(output, rule)

    @test output[1] == init
    @test output[2] ≈ [10.1     40.0029  70.0006;
                       20.02    50.0015  80.0004;
                       30.0067  60.0009  90.0003] atol=1e-4
    @test output[3] ≈ [10.198   40.0059  70.0011;
                       20.0399  50.0031  80.0008;
                       30.0133  60.0018  90.0005] atol=1e-4
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

    @test output[1] == LC50init
    @test output[2] ≈ [10.1923  40.0305  70.0097;
                       20.05    50.0055  80.0015;
                       30.0     60.0     90.0] atol=1e-4
    @test output[3] ≈ [10.3807  40.0609  70.0193
                       20.0998  50.011   80.0029
                       30.0     60.0     90.0] atol=1e-4
end

@testset "selection gradient allele frequency and phenotype on LC50" begin

    pFreqinit =  [0.1 0.5;
                  0.5 0.8]

    LC50init =  [50.0 50.0;
                 5.0 5.0]

    LandeInit = [LandeVariable(pFreqinit[i,j], LC50init[i,j]) for i in 1:2, j in 1:2]

    output = ArrayOutput(LandeInit; tspan=1:3)
    rule = Ruleset(SelectionGradient1locusSurv(exposure = 10, hillcoefficient=2.0,  deviationPhenotype = 10.0, dominanceDegree = -1.0))
    sim!(output, rule)

    @test output[1] == LandeInit
    # @test output[2] ≈ [LandeVariable(0.101274, 50.023475) LandeVariable(1.3, 26.12);
    #                    LandeVariable(0.503846, 50.0770414) LandeVariable(1.34272, 20.187476)] atol=1e-4
        
    # @test output[2] ≈ reshape([ # Have to use reshape to construct this
    #                         LandeVariable(0.101274, 50.023475), LandeVariable(1.3, 26.12),
    #                         LandeVariable(0.503846, 50.0770414), LandeVariable(1.34272, 20.187476)
    #                     ], (2, 2)) atol=1e-1

end

@testset "selection gradient allele frequency and phenotype on LC50 Map" begin

    pFreqinit =  [0.1 0.5;
                  0.5 0.8]

    LC50init =  [50.0 50.0;
                 5.0 5.0]

    LandeInit = [LandeVariable(pFreqinit[i,j], LC50init[i,j]) for i in 1:2, j in 1:2]

    exposure = [50.0 50.0;
                10.0 10.0]

    output = ArrayOutput(LandeInit; tspan=1:3, aux=(exposure=exposure,))
    rule = Ruleset(SelectionGradient1locusSurvMap(layerkey=Val(:exposure), hillcoefficient=2.0,  deviationPhenotype = 10.0, dominanceDegree = -1.0))
    sim!(output, rule)

    @test output[1] == LandeInit

end
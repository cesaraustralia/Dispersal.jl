using Dispersal, Test

@testset "selection gradient on LC50" begin
    init =  [10.0 40.0 70.0;
             20.0 50.0 80.0;
             30.0 60.0 90.0]
    output = ArrayOutput(init; tspan=1:3)
    rule = SelectionGradientSurv(; 
        exposure=10, hillcoefficient=2.0, additive_genetic_variance=1.0
    )
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
    lc50init =  [10.0 40.0 70.0;
                 20.0 50.0 80.0;
                 30.0 60.0 90.0]
    exposure = [50.0 50.0 50.0;
                20.0 20.0 20.0;
                0.0 0.0 0.0]

    output = ArrayOutput(lc50init; tspan=1:3, aux=(exposure=exposure,))
    rule = SelectionGradientSurvMap(
        exposurekey=Val(:exposure), hillcoefficient=2.0, additive_genetic_variance=1.0
    )
    sim!(output, rule)

    @test output[1] == lc50init
    @test output[2] ≈ [10.1923  40.0305  70.0097;
                       20.05    50.0055  80.0015;
                       30.0     60.0     90.0] atol=1e-4
    @test output[3] ≈ [10.3807  40.0609  70.0193
                       20.0998  50.011   80.0029
                       30.0     60.0     90.0] atol=1e-4
end

@testset "selection gradient allele frequency and phenotype on LC50" begin
    pfreqinit =  [0.1 0.5;
                  0.5 0.8]
    lc50init =  [50.0 50.0;
                 5.0 5.0]
    landeinit = [LandeVariable(pfreqinit[i,j], lc50init[i,j]) for i in 1:2, j in 1:2]

    output = ArrayOutput(landeinit; tspan=1:3)
    rule = SelectionGradient1locusSurv(
        exposure=10, hillcoefficient=2.0,  deviation_phenotype=10.0, dominance_degree=-1.0
    )
    sim!(output, rule)

    @test output[1] == landeinit
    # @test output[2] ≈ [LandeVariable(0.101274, 50.023475) LandeVariable(1.3, 26.12);
    #                    LandeVariable(0.503846, 50.0770414) LandeVariable(1.34272, 20.187476)] atol=1e-4
        
    # @test output[2] ≈ reshape([ # Have to use reshape to construct this
    #                         LandeVariable(0.101274, 50.023475), LandeVariable(1.3, 26.12),
    #                         LandeVariable(0.503846, 50.0770414), LandeVariable(1.34272, 20.187476)
    #                     ], (2, 2)) atol=1e-1
end

@testset "selection gradient allele frequency and phenotype on LC50 Map" begin
    pfreqinit = [0.1 0.5;
                 0.5 0.8]
    lc50init = [50.0 50.0;
                5.0 5.0]
    landeinit = [LandeVariable(pfreqinit[i,j], lc50init[i,j]) for i in 1:2, j in 1:2]
    exposure = [50.0 50.0;
                10.0 10.0]
    output = ArrayOutput(landeinit; tspan=1:3, aux=(exposure=exposure,))
    rule = SelectionGradient1locusSurvMap(
        exposurekey=Val(:exposure), hillcoefficient=2.0,  
        deviation_phenotype=10.0, dominance_degree=-1.0
    )
    sim!(output, rule)

    @test output[1] == landeinit

    # second implementation:
    landeinit = (
        pfreq = [0.1 0.5;
                 0.5 0.8],
        lc50 = [50.0 50.0;
                5.0 5.0],
    )
    output = ArrayOutput(landeinit; tspan=1:3, aux=(exposure=exposure,))

    exposure = [50.0 50.0;
                10.0 10.0]
    grids = Tuple{:pfreq,:lc50}
    rule = SelectionGradientMapTuple{grids,grids}(
        exposurekey=:exposure,
        hillcoefficient=2.0,  
        deviation_phenotype=10.0,
        dominance_degree=-1.0
    )
    sim!(output, rule)

    @test output[1] == landeinit
    @test output[2].pfreq ≈ [0.11656 0.55;
                             1.3 1.34272] atol=1e-4
    @test output[2].lc50 ≈ [50.3108 51.02;
                              26.12 20.1875] atol=1e-4
end
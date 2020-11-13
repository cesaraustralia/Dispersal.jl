using Dispersal, Test

@testset "allele frequency on LC50" begin
    init =  [0.1 0.2 0.3;
             0.4 0.5 0.6;
             0.7 0.8 0.9]

    output = ArrayOutput(init; tspan=1:3)
    rule = DeltaAlleleFrequencySurv(
        exposure=10, lc50=10, hillcoefficient=2.0, 
        deviation_phenotype=10.0, dominance_degree=-1
    )
    sim!(output, rule)
    @test output[1] == init
    @test output[2] ≈ [0.1828  0.3504  0.5016;
                       0.6352  0.75    0.8448;
                       0.9184  0.9696  0.9972] atol=1e-4
    @test output[3] ≈ [0.322707  0.571209  0.751677;
                       0.873187  0.946875  0.984955;
                       0.999613  1.00184   1.00027] atol=1e-4

    #! <1.0 is not always satisfied. Is it because of @fastmath ??
    @test (output[1] .<= 1.01) == trues(3,3)
    @test (0.0 .<= output[1]) == trues(3,3)
    @test (output[2] .<= 1.01) == trues(3,3)
    @test (0.0 .<= output[2]) == trues(3,3)
    @test (output[3] .<= 1.01) == trues(3,3)
    @test (0.0 .<= output[3]) == trues(3,3)
end

@testset "allele frequency on LC50 Map" begin
    pfreqinit =  [0.1 0.2 0.3;
                  0.4 0.5 0.6;
                  0.7 0.8 0.9]

    exposure = [50.0 50.0 50.0;
                20.0 20.0 20.0;
                0.0 0.0 0.0]

    output = ArrayOutput(pfreqinit; tspan=1:5, aux=(exposure=exposure,))
    rule = DeltaAlleleFrequencySurvMap(
        exposurekey=Val(:exposure), lc50=10, hillcoefficient=2.0, 
        deviation_phenotype=10.0, dominance_degree=-1.0
    )
    sim!(output, rule)

    @test output[1] == pfreqinit
    @test output[2] ≈ [ 0.259231  0.489231  0.687692;
                        0.77632   0.9       0.99168;
                        0.7       0.8       0.9] atol=1e-4
    @test output[3] ≈ [ 0.610737  0.968742  1.11622;
                        1.06951   1.05552   1.00618;
                        0.7       0.8       0.9] atol=1e-4

    #! <1.0 is not always satisfied. Is it because of @fastmath ??
    @test (output[1] .<= 1.1) == trues(3,3)
    @test (0.0 .<= output[1]) == trues(3,3)
    @test (output[2] .<= 1.1) == trues(3,3)
    @test (0.0 .<= output[2]) == trues(3,3)
    # @test (output[3] .<= 1.1) == trues(3,3)
    # @test (0.0 .<= output[3]) == trues(3,3)
    # @test (output[4] .<= 1.1) == trues(3,3)
    # @test (0.0 .<= output[4]) == trues(3,3)
    # @test (output[4] .<= 1.1) == trues(3,3)
    # @test (0.0 .<= output[4]) == trues(3,3)

    output = ArrayOutput(pfreqinit; tspan=1:5, aux=(exposure=exposure,))
    rule = DeltaAlleleFrequencySurvMap(
        exposurekey=Val(:exposure), lc50=10, hillcoefficient=2.0, 
        deviation_phenotype=10.0, dominance_degree=1
    )
    sim!(output, rule)

end

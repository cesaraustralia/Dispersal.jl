using DynamicGrids, Dispersal, Test, Unitful, Dates
using DynamicGrids: applyrule, SimData, extent
using Unitful: d

@testset "loglogistic survival" begin
    init =  [10.0 40.0 70.0;
             20.0 50.0 80.0;
             30.0 60.0 90.0]

    test1 = [10.0 40.0 70.0;
             20.0 50.0 80.0;
             30.0 60.0 90.0]
    test2 = [ 3.0327   12.1308  21.2289;
              6.06539  15.1635  24.2616;
              9.09809  18.1962  27.2943]
    test3 = [ 0.919724  3.6789   6.43807;
              1.83945   4.59862  7.35779;
              2.75917   5.51835  8.27752]

    output = ArrayOutput(init; tspan=1:3)
    rule = SurvLogLogistic(exposure=100, lc50=50, hillcoefficient=1.2)
    sim!(output, rule)

    @test output[1] == test1
    @test output[2] ≈ test2 atol=1e-4
    @test output[3] ≈ test3 atol=1e-4

end

@testset "loglogistic survival with exposure from exposure map" begin
    init =  [10.0 40.0 70.0;
            20.0 50.0 80.0;
            30.0 60.0 90.0]

    exposure = [100.0 100.0 100.0;
                50.0  50.0  50.0;
                20.0  20.0  20.0]

    output = ArrayOutput(init; tspan=1:3, aux=(exposure=exposure,))
    rule = SurvLogLogisticMap(exposurekey=Val(:exposure))
    sim!(output, rule)

    @test output[1] == [10.0 40.0 70.0;
                        20.0 50.0 80.0;
                        30.0 60.0 90.0]
    @test output[2] ≈  [4.42688  17.7075  30.9882;
                        9.19701  22.9925  36.7881;
                        14.4803   28.9607  43.441] atol=1e-4
    @test output[3] ≈  [1.95973   7.83892  13.7181;
                        4.22925  10.5731   16.917;
                        6.98935  13.9787   20.968] atol=1e-4

    @test DynamicGrids.normalise.(output[3], 0.25, 4) ≈ [0.455928  1.0  1.0;
                                                          1.0       1.0  1.0;
                                                          1.0       1.0  1.0] atol=1e-4

    output = ArrayOutput(init; tspan=1:3, aux=(exposure=exposure,))
    rule = SurvLogLogisticMap(
        exposurekey=Val(:exposure), 
        lc50=10
    )
    sim!(output, rule)

    @test output[1] == [10.0 40.0 70.0;
                        20.0 50.0 80.0;
                        30.0 60.0 90.0]
    @test output[2] ≈ [4.42688  17.7075  30.9882;
                       9.19701  22.9925  36.7881;
                      14.4803   28.9607  43.441] atol=1e-4
    @test output[3] ≈ [1.95973   7.83892  13.7181
                       4.22925  10.5731   16.917
                       6.98935  13.9787   20.968] atol=1e-4
end

@testset "loglogistic survival with ArrayOutput layer" begin
    init =  [10.0 40.0 70.0;
            20.0 50.0 80.0;
            30.0 60.0 90.0]
    suit = [100.0 100.0 100.0;
            50.0  50.0  50.0;
            20.0  20.0  20.0]
    # suit = output
    # outputStack = ArrayOutput(init; tspan=1:3, aux=(suit=suit,index=1:3))
    # rule = SurvLogLogisticMap(layerkey=:suit, timeindex=:index)
    # sim!(outputStack, rule)
end

@testset "survival map masking" begin
    init = [1.0 1.0 1.0;
            1.0 1.0 1.0;
            1.0 1.0 1.0]

    exposure =  [1.0 1.0 2.0;
                 2.0 2.0 0.5;
                 1.0 1.0 0.5]

    test1 = [1.0 1.0 1.0;
             1.0 1.0 1.0;
             1.0 1.0 1.0]
    test2 = [0.0 0.0 1.0;
             1.0 1.0 0.0;
             0.0 0.0 0.0]

    output = ArrayOutput(init; tspan=1:3, aux=(exposure=exposure,))
    maskrule = MaskSurvMap(exposurekey=Val(:exposure), threshold=1.1)
    data = SimData(extent(output), Ruleset(maskrule))

    @test applyrule(data, maskrule, 5, (1, 1)) === 0
    @test applyrule(data, maskrule, 5, (2, 2)) === 5
    @test applyrule(data, maskrule, 5.0, (1, 1)) === 0.0
    @test applyrule(data, maskrule, 5.0, (2, 2)) === 5.0

    sim!(output, maskrule)

    @test output[1] == test1
    @test output[2] == test2
end

@testset "Test double layers ExactLogisticGrowthMap2" begin
    popsizeinit = [1.0 4.0 7.0;
                   2.0 5.0 8.0;
                   3.0 6.0 9.0]
    exposure = repeat([0.0  0.0  0.0;
                       20.0 20.0 20.0;
                       50.0 50.0 50.0], inner=(1, 1, 3))
    lc50 = cat([0.0 0.0 0.0;
                0.0 0.0 0.0;
                0.0 0.0 0.0],
               [20.0 20.0 20.0;
                20.0 20.0 20.0;
                20.0 20.0 20.0],
               [10.0 10.0 10.0;
                10.0 10.0 10.0;
                10.0 10.0 10.0]; dims=3)

    popsizegrids = ArrayOutput(popsizeinit; tspan=1:6, 
        aux=(exposure=exposure, lc50=lc50)
    )
    survrule = SurvLogLogisticMap2(exposurekey=:exposure, lc50key=:lc50);
    sim!(popsizegrids, survrule);
end

@testset "Test double grids GrowthSurvLogLogisticMap3" begin
    init = (
        pop1=[0.5 0.0 0.0;
              0.0 0.0 1.0;
              0.0 0.0 0.0],
        pop2=[0.5 0.0 0.0;
              0.0 1.0 0.0;
              0.0 0.0 0.0],
        pop3=[0.5 0.0 0.0;
              0.0 1.0 0.0;
              0.0 1.0 0.0],
    )
    exposure = repeat([ 0.0  0.0  0.0;
                       20.0 20.0 20.0;
                       50.0 50.0 50.0], inner=(1, 1,3))

    output = ArrayOutput(init; tspan=1:3, aux=(exposure=exposure,))
    grids = Tuple{:pop1,:pop2,:pop3}
    rule = GrowthSurvLogLogisticMap3{grids,grids}(
        exposurekey=:exposure,
        intrinsicrate=1.5,
        hillcoefficient=2.5,
        lc50=(50, 40, 30),
        carrycap = 5.0,
    )
    GrowthSurvLogLogisticMap3{grids,grids}(1,2,3,4,4,5,6,7)

    sim!(output, rule)

    @test output[1] == init
    @test output[2].pop1 ≈ [.525 0. 0.;
                            0. 0. 1.08973;
                            0. 0. 0.] atol=1e-4
    @test output[2].pop2 ≈ [.525 0. 0.;
                            0. 0.764801  0.;
                            0. 0. 0.] atol=1e-4
    @test output[2].pop3 ≈ [.525 0. 0.;
                            0. 0.660363 0.;
                            0. 0.26166 0.] atol=1e-4

    # ----------------- SIMPLER SAME FUNCTION
    ruleSurv = SurvLogLogisticMapTuple{grids,grids}(
        exposurekey=:exposure,
        hillcoefficient=2.5,
        lc50=(50, 40, 30),
    )
    ruleGrowth = GrowthMapTuple{grids,grids}(
        carrycap=5.0,
        intrinsicrate=1.5,
    )
    rule = Ruleset(ruleGrowth, ruleSurv)

    sim!(output, rule)

    @test output[1] == init
    @test output[2].pop1 ≈ [.525 0. 0.;
                            0. 0. 1.08973;
                            0. 0. 0.] atol=1e-4
    @test output[2].pop2 ≈ [.525 0. 0.;
                            0. 0.764801  0.;
                            0. 0. 0.] atol=1e-4
    @test output[2].pop3 ≈ [.525 0. 0.;
                            0. 0.660363 0.;
                            0. 0.26166 0.] atol=1e-4
end

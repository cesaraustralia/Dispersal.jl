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
    rule = Ruleset(SurvLogLogistic(exposure = 100, LC50=50, hillcoefficient = 1.2))
    sim!(output, rule)

    @test output[1] == test1
    @test output[2] ≈ test2 atol=1e-4
    @test output[3] ≈ test3 atol=1e-4

end

@testset "loglogistic survival with exposure from exposure map" begin
    init =  [10.0 40.0 70.0;
            20.0 50.0 80.0;
            30.0 60.0 90.0]

    suit = [100.0 100.0 100.0;
            50.0  50.0  50.0;
            20.0  20.0  20.0]

    output = ArrayOutput(init; tspan=1:3, aux=(suit=suit,))
    output.extent
    rule = Ruleset(SurvLogLogisticMap(layerkey=Val(:suit)))
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

    output = ArrayOutput(init; tspan=1:3, aux=(suit=suit,))
    rule = Ruleset(SurvLogLogisticMap(layerkey=Val(:suit), LC50=10))
    sim!(output, rule)

    @test output[1] == test1
    @test output[2] ≈ [4.42688  17.7075  30.9882;
                       9.19701  22.9925  36.7881;
                      14.4803   28.9607  43.441] atol=1e-4
    @test output[3] ≈ [1.95973   7.83892  13.7181
                       4.22925  10.5731   16.917
                       6.98935  13.9787   20.968] atol=1e-4
end

@testset "survival map masking" begin
    init = [1.0 1.0 1.0;
        1.0 1.0 1.0;
        1.0 1.0 1.0]

    suit =  [1.0 1.0 2.0;
        2.0 2.0 0.5;
        1.0 1.0 0.5]

    test1 = [1.0 1.0 1.0;
        1.0 1.0 1.0;
        1.0 1.0 1.0]

    test2 = [0.0 0.0 1.0;
        1.0 1.0 0.0;
        0.0 0.0 0.0]

    output = ArrayOutput(init; tspan=1:3, aux=(suit=suit,))
    maskrule = MaskSurvMap(layerkey=Val(:suit), threshold=1.1)
    ruleset = Ruleset(maskrule)
    data = SimData(extent(output), ruleset)

    @test applyrule(data, maskrule, 5, (1, 1)) === 0
    @test applyrule(data, maskrule, 5, (2, 2)) === 5

    @test applyrule(data, maskrule, 5.0, (1, 1)) === 0.0
    @test applyrule(data, maskrule, 5.0, (2, 2)) === 5.0

    sim!(output, Ruleset(maskrule))

    @test output[1] == test1
    @test output[2] == test2
end

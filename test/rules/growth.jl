using DynamicGrids, Dispersal, Test, Unitful, Dates
using DynamicGrids: applyrule, SimData, extent
using Unitful: d

@testset "exponential growth" begin
    init =  [1.0 4.0 7.0;
             2.0 5.0 8.0;
             3.0 6.0 9.0]

    test1 = [ 1.0  4.0  7.0;
              2.0  5.0  8.0;
              3.0  6.0  9.0]
    test2 = [ 2.0  8.0 14.0;
              4.0 10.0 16.0;
              6.0 12.0 18.0]
    test3 = [ 4.0 16.0 28.0;
              8.0 20.0 32.0;
             12.0 24.0 36.0]

    output = ArrayOutput(init; tspan=1:3)
    rule = Ruleset(ExponentialGrowth(intrinsicrate=log(2.0)))
    sim!(output, rule)

    @test output[1] == test1
    @test output[2] == test2
    @test output[3] == test3

    output = ArrayOutput(init; tspan=1d:1d:3d)
    rule = Ruleset(ExponentialGrowth(intrinsicrate=log(2.0), timestep=1d); timestep=1d)
    sim!(output, rule)

    @test output[1] == test1
    @test output[2] == test2
    @test output[3] == test3

    output = ArrayOutput(init; tspan=DateTime(2001,1,1):Day(5):DateTime(2001,1,15))
    rule = Ruleset(ExponentialGrowth(intrinsicrate=log(2.0)/5, timestep=Day(1)); timestep=Day(5))
    sim!(output, rule)

    @test output[1] == test1
    @test output[2] == test2
    @test output[3] == test3
end

@testset "exponential growth with rate from growth map" begin
    init = [1.0 1.0 1.0;
            1.0 1.0 1.0;
            1.0 1.0 1.0]

    suit =  log.([1.0 1.0 2.0;
                  2.0 1.0 0.5;
                  1.0 1.0 0.5])

    output = ArrayOutput(init; tspan=1:3, aux=(suit=suit,))
    output.extent
    rule = Ruleset(ExponentialGrowthMap(auxkey=Val(:suit)))
    sim!(output, rule)

    @test output[1] == [1.0 1.0 1.0;
                        1.0 1.0 1.0;
                        1.0 1.0 1.0]
    @test output[2] ≈  [1.0 1.0 2.0;
                        2.0 1.0 0.5;
                        1.0 1.0 0.5]
    @test output[3] ≈  [1.0 1.0 4.0;
                        4.0 1.0 0.25;
                        1.0 1.0 0.25]

    @test DynamicGrids.normalise.(output[3], 0.25, 4) == [0.2  0.2  1.0;
                                                          1.0  0.2  0.0;
                                                          0.2  0.2  0.0]

end

@testset "logistic growth" begin

    init =  [1.0 4.0 7.0;
             2.0 5.0 8.0;
             3.0 6.0 9.0]

    test1 = [1.0  4.0  7.0;
             2.0  5.0  8.0;
             3.0  6.0  9.0]

    test2 = [1.81818  5.71429  8.23529;
             3.33333  6.66667  8.88889;
             4.61538  7.5      9.47368]

    test3 = [3.07692  7.27273  9.03226;
             5.0      8.0      9.41176;
             6.31579  8.57143  9.72973]
    rule = Ruleset(LogisticGrowth(intrinsicrate=log(2.0), carrycap=10))
    output = ArrayOutput(init; tspan=1:3)

    sim!(output, rule)

    @test output[1] == test1
    @test output[2] ≈ test2 atol=1e-4
    @test output[3] ≈ test3 atol=1e-4
end

@testset "logistic growth with rate from a growth map" begin

    init = [1.0 1.0 1.0;
            1.0 1.0 1.0;
            1.0 1.0 1.0]

    test1 = [1.0 1.0 1.0;
             1.0 1.0 1.0;
             1.0 1.0 1.0]

    test2 = [1.0     1.0  1.81818;
             1.81818 1.0  0.5    ;
             1.0     1.0  0.5    ]

    test3 = [1.0      1.0  3.07692;
             3.07692  1.0  0.25   ;
             1.0      1.0  0.25   ]

    suit =  log.([1.0 1.0 2.0;
                  2.0 1.0 0.5;
                  1.0 1.0 0.5])

    output = ArrayOutput(init; tspan=1:3, aux=(suit=suit,))
    rule = Ruleset(LogisticGrowthMap(auxkey=Val(:suit), carrycap=10))
    sim!(output, rule)

    @test output[1] == test1
    @test output[2] ≈ test2 atol=1e-4
    @test output[3] ≈ test3 atol=1e-4

end

@testset "growth map masking" begin
    init =  [1.0 1.0 1.0;
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
    maskrule = MaskGrowthMap(auxkey=Val(:suit), threshold=1.1)
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

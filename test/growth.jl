using Cellular, Dispersal, Test
using Cellular: rule, FrameData

data = FrameData(nothing, nothing, nothing, nothing, nothing, 1)

@testset "exponential growth" begin
    init =  [1.0 4.0 7.0;
             2.0 5.0 8.0;
             3.0 6.0 9.0]

    output = ArrayOutput(init, 3)
    model = Models(ExactExponentialGrowth(intrinsicrate = log(2.0)))
    sim!(output, model, init; tstop=3)

    @test output[1] == [ 1.0  4.0  7.0;
                         2.0  5.0  8.0;
                         3.0  6.0  9.0]
    @test output[2] == [ 2.0  8.0 14.0;
                         4.0 10.0 16.0;
                         6.0 12.0 18.0]
    @test output[3] == [ 4.0 16.0 28.0;
                         8.0 20.0 32.0;
                        12.0 24.0 36.0]
end

@testset "exponential growth with rate from suitability layer" begin
    init = [1.0 1.0 1.0;
            1.0 1.0 1.0;
            1.0 1.0 1.0]

    suit =  log.([1.0 1.0 2.0;
                  2.0 1.0 0.5;
                  1.0 1.0 0.5])

    output = ArrayOutput(init, 3)
    model = Models(SuitabilityExactExponentialGrowth(layers = suit))
    sim!(output, model, init; tstop=3)

    @test output[1] == [1.0 1.0 1.0;
                        1.0 1.0 1.0;
                        1.0 1.0 1.0]
    @test output[2] ≈  [1.0 1.0 2.0;
                        2.0 1.0 0.5;
                        1.0 1.0 0.5]
    @test output[3] ≈  [1.0 1.0 4.0;
                        4.0 1.0 0.25;
                        1.0 1.0 0.25]

    @test Cellular.normalize_frame(output[3], 0.25, 4) == [0.2  0.2  1.0;
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

    model = Models(ExactLogisticGrowth(intrinsicrate = log(2.0), carrycap = 10))
    output = ArrayOutput(init, 3)

    sim!(output, model, init; tstop=3)

    @test output[1] == test1
    @test output[2] ≈ test2 atol=1e-4
    @test output[3] ≈ test3 atol=1e-4
end

@testset "logistic growth with rate from suitability layer" begin

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

    output = ArrayOutput(init, 3)
    model = Models(SuitabilityExactLogisticGrowth(layers = suit, carrycap = 10))
    sim!(output, model, init; tstop=3)

    @test output[1] == test1
    @test output[2] ≈ test2 atol=1e-4
    @test output[3] ≈ test3 atol=1e-4

end

@testset "suitability masking" begin
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

    output = ArrayOutput(init, 3)
    mask = SuitabilityMask(layers=suit, threshold=1.1)

    @test rule(mask, data, 5, (1, 1)) === 0
    @test rule(mask, data, 5, (2, 2)) === 5

    @test rule(mask, data, 5.0, (1, 1)) === 0.0
    @test rule(mask, data, 5.0, (2, 2)) === 5.0

    sim!(output, Models(mask), init; tstop=3)

    @test output[1] == test1
    @test output[2] == test2
end


@testset "growth models" begin

    global init =  setup([1.0 4.0 7.0;
                          2.0 5.0 8.0;
                          3.0 6.0 9.0])

    global output = ArrayOutput(init, 3)


    @testset "exponential growth" begin
        global model = Models(ExactExponentialGrowth(intrinsicrate = log(2.0), timestep = 1))
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

        minmaxrange = 0.0, 40.0

        global init = setup([1.0 1.0 1.0;
                             1.0 1.0 1.0;
                             1.0 1.0 1.0])
        global output = ArrayOutput(init, 3)

        global suit =  log.(setup([1.0 1.0 2.0;
                                   2.0 1.0 0.5;
                                   1.0 1.0 0.5]))
        global model = Models(SuitabilityExactExponentialGrowth(layers = suit))

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

end

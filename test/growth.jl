
@testset "growth models" begin

    global init =  setup([1.0 4.0 7.0;
                          2.0 5.0 8.0;
                          3.0 6.0 9.0])
    
    global output = ArrayOutput(init, 3)


    @testset "exponential growth" begin
        global model = Models(ExponentialGrowth(2.0))
        sim!(output, model, init; tstop=3)
        @test output[1] == [1.0 4.0 7.0;
                            2.0 5.0 8.0;
                            3.0 6.0 9.0]
        @test output[2] == [2.0  8.0 14.0;
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

        global suit =  setup([1.0 4.0 7.0;
                              2.0 5.0 0.5;
                              3.0 6.0 -1.0])
        global model = Models(SuitabilityExponentialGrowth(suit, minmaxrange...))

        sim!(output, model, init; tstop=3)
        @test output[1] == [1.0 1.0 1.0;
                            1.0 1.0 1.0;
                            1.0 1.0 1.0]
        @test output[2] ≈  [1.0 4.0 7.0;
                            2.0 5.0 0.5;
                            3.0 6.0 0.0]
        @test output[3] ≈  [1.0 16.0 40.0;
                            4.0 25.0 0.25;
                            9.0 36.0 0.0]

        @test Cellular.normalize_frame(output[3], minmaxrange...) == [0.025 0.4 1.0;
                                                                      0.1 0.625 0.00625;
                                                                      0.225 0.9 0.0] 

    end

end

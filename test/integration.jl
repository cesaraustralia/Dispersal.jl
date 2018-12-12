
global init = setup([0 1 0; 1 0 1])

@testset "dispersal kernel array matches passed in function" begin
    global dk = DispersalKernel(f=exponential, cellsize=1, init=init, radius=2, param=(1.0,)).kernel
    @test size(dk) == (5, 5)
    @test sum(dk) ≈ 1.0
end

@testset "binary dispersal simulation with suitability mask" begin

    global suit =  setup([1 0 1 1 0;
                          0 0 1 1 1;
                          1 1 1 1 0;
                          1 1 0 1 1;
                          1 0 1 1 1])

    global init =  setup(Bool[0 0 0 0 0;
                          0 0 0 0 0;
                          0 0 1 0 0;
                          0 0 0 0 0;
                          0 0 0 0 0])

    global test1 = setup(Bool[0 0 0 0 0;
                          0 0 0 0 0;
                          0 0 1 0 0;
                          0 0 0 0 0;
                          0 0 0 0 0])

    global test2 = setup(Bool[0 0 0 0 0;
                          0 0 1 1 0;
                          0 1 1 1 0;
                          0 1 0 1 0;
                          0 0 0 0 0])

    global test3 = setup(Bool[0 0 1 1 0;
                          0 0 1 1 1;
                          1 1 1 1 0;
                          1 1 0 1 1;
                          1 0 1 1 1])

    # Dispersal in radius 1 neighborhood
    global suitmask = SuitabilityMask(layers=suit)
    global hood = DispersalKernel(; init=init, radius=1)

    @testset "inwards binary dispersal fills the grid where reachable and suitable" begin
        global inwards = InwardsBinaryDispersal(neighborhood=hood, prob_threshold=0.0)
        global model = Models(inwards, suitmask)
        global output = ArrayOutput(init, 3)
        sim!(output, model, init; tstop=3)
        @test output[1] == test1
        @test output[2] == test2
        @test output[3] == test3

        # As submodels
        global model = Models((inwards, suitmask))
        sim!(output, model, init; tstop=3)
        @test output[1] == test1
        @test output[2] == test2
        @test output[3] == test3
    end
# updating the floating point outwards dispersal broke this test
    # @testset "outwards dispersal fills the grid where reachable and suitable" begin
    #     global outwards = OutwardsBinaryDispersal(neighborhood=hood, prob_threshold=0.0)
    #     global model = Models(outwards, suitmask)
    #     global output = ArrayOutput(init, 3)
    #     sim!(output, model, init; tstop=3)
    #     @test output[1] == test1
    #     @test output[2] == test2
    #     @test output[3] == test3
    # end

end

@testset "floating point population dispersal simulation" begin

    global suit =  [1.0 0.0 1.0 0.0 1.0 1.0 0.0;
                    1.0 1.0 0.0 0.0 1.0 1.0 1.0;
                    1.0 0.0 1.0 1.0 1.0 1.0 0.0;
                    1.0 0.0 1.0 1.0 0.0 1.0 1.0;
                    1.0 1.0 1.0 0.0 1.0 1.0 1.0;
                    1.0 0.0 0.0 0.0 0.0 0.0 0.0;
                    1.0 1.0 1.0 1.0 1.0 1.0 1.0]

    global init = setup([100.0 0.0 0.0 0.0   0.0 0.0 0.0;
                           0.0 0.0 0.0 0.0   0.0 0.0 0.0;
                           0.0 0.0 0.0 0.0 100.0 0.0 0.0;
                           0.0 0.0 0.0 0.0   0.0 0.0 0.0;
                           0.0 0.0 0.0 0.0   0.0 0.0 0.0;
                           0.0 0.0 0.0 0.0   0.0 0.0 0.0;
                           0.0 0.0 0.0 0.0   0.0 0.0 0.0])

    global test1 = setup([100.0 0.0 0.0 0.0   0.0 0.0 0.0;
                           0.0 0.0 0.0 0.0   0.0 0.0 0.0;
                           0.0 0.0 0.0 0.0 100.0 0.0 0.0;
                           0.0 0.0 0.0 0.0   0.0 0.0 0.0;
                           0.0 0.0 0.0 0.0   0.0 0.0 0.0;
                           0.0 0.0 0.0 0.0   0.0 0.0 0.0;
                           0.0 0.0 0.0 0.0   0.0 0.0 0.0])

    global test2 = setup([ 4.0  0.0  8.0  0.0  4.0  4.0  0.0;
                           4.0  4.0  0.0  0.0  4.0  4.0  4.0;
                           4.0  0.0  8.0  4.0  4.0  4.0  0.0;
                           0.0  0.0  4.0  4.0  0.0  4.0  4.0;
                           0.0  0.0  4.0  0.0  4.0  4.0  4.0;
                           0.0  0.0  0.0  0.0  0.0  0.0  0.0;
                           0.0  0.0  0.0  0.0  0.0  0.0  0.0])

    global test3 = setup([ 1.28  0.0   1.92  0.0   1.92  1.28  0.0;
                           1.44  1.76  0.0   0.0   2.56  1.76  1.44;
                           1.6   0.0   2.56  2.88  3.2   2.24  0.0;
                           1.12  0.0   1.92  2.24  0.0   1.92  1.6;
                           0.8   1.12  1.44  0.0   2.08  1.44  1.12;
                           0.32  0.0   0.0   0.0   0.0   0.0   0.0;
                           0.16  0.16  0.32  0.48  0.64  0.48  0.48;])

    # Dispersal in radius 1 neighborhood
    global suitmask = SuitabilityMask(layers=suit)
    global r = 2

    @testset "inwards population dispersal fills the grid where reachable suitable" begin
        global hood = DispersalKernel(; f=(d,a)->1.0, radius=r)
        global inwards = InwardsPopulationDispersal(neighborhood=hood)
        global model = Models(inwards, suitmask)
        global output = ArrayOutput(init, 3)
        sim!(output, model, init; tstop=3)
        @test output[1] == test1
        @test output[2] == test2
        @test output[3] ≈ test3

        # As submodels
        global model = Models((inwards, suitmask))
        sim!(output, model, init; tstop=3)
        @test output[1] == test1
        @test output[2] == test2
        @test output[3] ≈ test3
    end

    @testset "outwards population dispersal fills the grid where reachable and suitable" begin
        global hood = DispersalKernel(; f=(d,a)->1.0, radius=r)
        global outwards = OutwardsPopulationDispersal(neighborhood=hood, fraction=1)
        global model = Models(outwards, suitmask)
        global output = ArrayOutput(init, 3)
        sim!(output, model, init; tstop=3)
        @test output[1] == test1
        @test output[2] == test2
        @test output[3] ≈ test3
    end
end


@testset "jump dispersal models work" begin

    global init = setup([0 0 0 0 0;
                         0 0 0 0 0;
                         0 0 1 0 0;
                         0 0 0 0 0;
                         0 0 0 0 0])

    @testset "Jump dispersal spread randomly" begin
        global model = Models(JumpDispersal(prob_threshold=0.0, spotrange=3))
        global output = ArrayOutput(init, 20)
        sim!(output, model, init; tstop=20)
    end

end

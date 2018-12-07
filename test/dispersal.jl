
global init = setup([0 1 0; 1 0 1])

@testset "dispersal kernel array matches passed in function" begin
    global dk = DispersalNeighborhood(f=exponential, cellsize=1, init=init, radius=2, param=(1.0,)).kernel
    @test size(dk) == (5, 5)
    @test sum(dk) ≈ 1.0
end

@testset "dispersal nieghborhood sum matches the passed-in kernel function" begin
    global hood = DispersalNeighborhood(f=exponential, init=init, cellsize=1, radius=2)
    global state = 0
    global t = 0

    global source = setup([1 0 0 0 0;
                           0 0 0 1 0;
                           0 1 0 0 1;
                           0 0 0 1 0;
                           0 1 0 0 0])

    @testset "neighborhood sum matches grid * kernel sum for same-sized grid" begin
        modeldata = Cellular.prepare((InwardsBinaryDispersal(),), init)
        data = Cellular.FrameData(source, deepcopy(source), 1, 1, modeldata)
        # global cc = neighbors(hood, nothing, data, state)
        # @test cc ≈ sum(source .* hood.kernel)
    end

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
    global layers = SuitabilityLayer(suit)
    global suitmask = SuitabilityMask()
    global hood = DispersalNeighborhood(; init=init, radius=1)

    @testset "inwards binary dispersal fills the grid where reachable and suitable" begin
        global inwards = InwardsBinaryDispersal(neighborhood=hood, prob_threshold=0.0)
        global model = Models(inwards, suitmask)
        global output = ArrayOutput(init)
        sim!(output, model, init, layers; tstop=3)
        @test output[1] == test1
        @test output[2] == test2
        @test output[3] == test3
    end

    @testset "outwards dispersal fills the grid where reachable and suitable" begin
        global outwards = OutwardsBinaryDispersal(neighborhood=hood, prob_threshold=0.0)
        global model = Models(outwards, suitmask)
        global output = ArrayOutput(init)
        sim!(output, model, init, layers; tstop=3)
        @test output[1] == test1
        @test output[2] == test2
        @test output[3] == test3
    end

end

@testset "floating point population dispersal simulation" begin

    setup(x) = ScalableMatrix(x, 0.0, 100.0)

    global suit =  [1.0 0.0 1.0 0.0 1.0 1.0 0.0;
                    1.0 1.0 0.0 0.0 1.0 1.0 1.0;
                    1.0 0.0 1.0 1.0 1.0 1.0 0.0;
                    1.0 0.0 1.0 1.0 0.0 1.0 1.0;
                    1.0 1.0 1.0 0.0 1.0 1.0 1.0;
                    1.0 0.0 0.0 0.0 0.0 0.0 0.0;
                    1.0 1.0 1.0 1.0 1.0 1.0 1.0]

    global init =  setup([1.0 0.0 0.0 0.0 0.0 0.0 0.0;
                          0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                          0.0 0.0 0.0 0.0 1.0 0.0 0.0;
                          0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                          0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                          0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                          0.0 0.0 0.0 0.0 0.0 0.0 0.0])

    global test1 = setup([1.0 0.0 0.0 0.0 0.0 0.0 0.0;
                          0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                          0.0 0.0 0.0 0.0 1.0 0.0 0.0;
                          0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                          0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                          0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                          0.0 0.0 0.0 0.0 0.0 0.0 0.0])

    global test2 = setup([2.0 0.0 2.0 0.0 1.0 1.0 0.0;
                          1.0 1.0 0.0 0.0 1.0 1.0 1.0;
                          1.0 0.0 2.0 1.0 2.0 1.0 0.0;
                          0.0 0.0 1.0 1.0 0.0 1.0 1.0;
                          0.0 0.0 1.0 0.0 1.0 1.0 1.0;
                          0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                          0.0 0.0 0.0 0.0 0.0 0.0 0.0])

    global test3 = setup([11.0  0.0 16.0  0.0 14.0  10.0  0.0;
                          11.0 13.0  0.0  0.0 18.0  13.0 11.0;
                          12.0  0.0 20.0 20.0 23.0  16.0  0.0;
                          7.0   0.0 14.0 16.0  0.0  14.0 12.0;
                          5.0   7.0 11.0  0.0 15.0  11.0  9.0;
                          2.0   0.0  0.0  0.0  0.0   0.0  0.0;
                          1.0   1.0  2.0  3.0  4.0   3.0  3.0])

    # Dispersal in radius 1 neighborhood
    global layers = SuitabilityLayer(suit)
    global suitmask = SuitabilityMask()
    global r = 2

    @testset "inwards population dispersal fills the grid where reachable suitable" begin
        global hood = DispersalNeighborhood(; f=(d,a)->1.0, radius=r)
        global inwards = InwardsPopulationDispersal(neighborhood=hood, fraction=(2r+1)^2)
        global model = Models(inwards, suitmask)
        global output = ArrayOutput(init)
        sim!(output, model, init, layers; tstop=3)
        @test output[1] == test1
        @test output[2] == test2
        @test output[3] ≈ test3
    end

    @testset "outwards population dispersal fills the grid where reachable and suitable" begin
        global hood = DispersalNeighborhood(; f=(d,a)->1.0, radius=r)
        global outwards = OutwardsPopulationDispersal(neighborhood=hood, fraction=(2r+1)^2)
        global model = Models(outwards, suitmask)
        global output = ArrayOutput(init)
        sim!(output, model, init, layers; tstop=3)
        @test output[1] == test1
        @test output[2] == test2
        @test output[3] == test3
    end
end


@testset "jump dispersal models work" begin

    global init = setup([0 0 0 0 0;
                         0 0 0 0 0;
                         0 0 1 0 0;
                         0 0 0 0 0;
                         0 0 0 0 0])

    @testset "Jump dispersal spread randomly" begin
        global layers = SuitabilityLayer(suit)
        global model = Models(JumpDispersal(prob_threshold=0.0, spotrange=3))
        global output = ArrayOutput(init)
        sim!(output, model, init, layers; tstop=20)
    end

end

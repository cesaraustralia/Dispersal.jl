
global init = setup([0 1 0; 1 0 1])

@testset "dispersal kernel array matches passed in function" begin
    global dk = DispersalNeighborhood(dir=:inwards, f=exponential, cellsize=1, init=init, radius=2, param=(1.0,)).kernel
    @test size(dk) == (5, 5)
    @test sum(dk) ≈ 1.0
end

@testset "dispersal nieghborhood sum matches the passed-in kernel function" begin
    global hood = DispersalNeighborhood(dir=:inwards, f=exponential, init=init, cellsize=1, radius=2)
    global state = 0
    global t = 0

    global source = setup([1 0 0 0 0;
                           0 0 0 1 0;
                           0 1 0 0 1;
                           0 0 0 1 0;
                           0 1 0 0 0])

    @testset "neighborhood sum matches grid * kernel sum for same-sized grid" begin
        data = Cellular.ModelData(1, source, deepcopy(source), 1)
        global cc = neighbors(hood, nothing, data, state, (3, 3))
        @test cc ≈ sum(source .* hood.kernel)
    end

end 


@testset "binary dispersal simulation with suitability mask" begin

    global suit =  setup([1 0 1 1 0;
                          0 0 1 1 1;
                          1 1 1 1 0;
                          1 1 0 1 1;
                          1 0 1 1 1])

    global init =  setup([0 0 0 0 0;
                          0 0 0 0 0;
                          0 0 1 0 0;
                          0 0 0 0 0;
                          0 0 0 0 0])

    global test1 = setup([0 0 0 0 0;
                          0 0 0 0 0;
                          0 0 1 0 0;
                          0 0 0 0 0;
                          0 0 0 0 0])

    global test2 = setup([0 0 0 0 0;
                          0 0 1 1 0;
                          0 1 1 1 0;
                          0 1 0 1 0;
                          0 0 0 0 0])

    global test3 = setup([0 0 1 1 0;
                          0 0 1 1 1;
                          1 1 1 1 0;
                          1 1 0 1 1;
                          1 0 1 1 1])

    # Dispersal in radius 1 neighborhood
    global layers = SuitabilityLayer(suit)
    global suitmask = SuitabilityMask()

    @testset "inwards binary dispersal fills the grid where reachable and suitable" begin
        global hood = DispersalNeighborhood(; dir=:inwards, init=init, radius=1)
        global inwards = InwardsBinaryDispersal(neighborhood=hood, prob_threshold=0.0)
        global model = Models(inwards, suitmask)
        global output = ArrayOutput(init)
        sim!(output, model, init, layers; time=3)
        @test output[1] == test1
        @test output[2] == test2
        @test output[3] == test3
    end

    @testset "outwards dispersal fills the grid where reachable and suitable" begin
        global hood = DispersalNeighborhood(; dir=:outwards,  radius=1)
        global outwards = OutwardsBinaryDispersal(neighborhood=hood, prob_threshold=0.0)
        global model = Models(outwards, suitmask)
        global output = ArrayOutput(init)
        sim!(output, model, init, layers; time=3)
        @test output[1] == test1
        @test output[2] == test2
        @test output[3] == test3
    end
end

@testset "floating point population dispersal simulation" begin

    global suit =  setup([1.0 0.0 1.0 1.0 0.0;
                          0.0 0.0 1.0 1.0 1.0;
                          1.0 1.0 1.0 1.0 0.0;
                          1.0 1.0 0.0 1.0 1.0;
                          1.0 0.0 1.0 1.0 1.0])

    global init =  setup([0.0 0.0 0.0 0.0 0.0;
                          0.0 0.0 0.0 0.0 0.0;
                          0.0 0.0 1.0 0.0 0.0;
                          0.0 0.0 0.0 0.0 0.0;
                          0.0 0.0 0.0 0.0 0.0])

    global test1 = setup([0.0 0.0 0.0 0.0 0.0;
                          0.0 0.0 0.0 0.0 0.0;
                          0.0 0.0 1.0 0.0 0.0;
                          0.0 0.0 0.0 0.0 0.0;
                          0.0 0.0 0.0 0.0 0.0])

    global test2 = setup([0.0 0.0 0.0 0.0 0.0;
                          0.0 0.0 1.0 1.0 0.0;
                          0.0 1.0 2.0 1.0 0.0;
                          0.0 1.0 0.0 1.0 0.0;
                          0.0 0.0 0.0 0.0 0.0])

    global test3 = setup([0.0 0.0 2.0  2.0 0.0;
                          0.0 0.0 7.0  6.0 2.0;
                          2.0 6.0 10.0 7.0 0.0;
                          2.0 5.0 0.0  5.0 2.0;
                          1.0 0.0 2.0  1.0 1.0])

    # Dispersal in radius 1 neighborhood
    global layers = SuitabilityLayer(suit)

    @testset "inwards population dispersal fills the grid where reachable suitable" begin
        global hood = DispersalNeighborhood(; dir=:inwards, f=(d,a)->1.0, radius=1)
        global inwards = InwardsPopulationDispersal(neighborhood=hood, fraction=9/1)
        global model = Models(inwards, suitmask)
        global output = ArrayOutput(init)
        sim!(output, model, init, layers; time=3)
        @test output[1] == test1
        @test output[2] == test2
        @test output[3] ≈ test3
    end

    @testset "outwards population dispersal fills the grid where reachable and suitable" begin
        global hood = DispersalNeighborhood(; dir=:outwards, f=(d,a)->1.0, radius=1)
        global outwards = OutwardsPopulationDispersal(neighborhood=hood, fraction=9/1)
        global model = Models(outwards, suitmask)
        global output = ArrayOutput(init)
        sim!(output, model, init, layers; time=3)
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
        sim!(output, model, init, layers; time=20)
    end

end

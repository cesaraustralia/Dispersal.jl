using Revise,
      Dispersal,
      Cellular,
      Test
using Dispersal: suitability, cyclic, sequence_interpolate, neighbors

setup(x) = x

# For manual testing on CUDA
# using CuArrays, CUDAnative
# setup(x) = CuArray(x)

@testset "suitability is 1.0 by default" begin
    @test suitability(nothing, (1, 1), 10) == 1.0
end

@testset "single layer suitability just returns the layer value" begin
    global suitlayer = SuitabilityLayer(setup([0.1 0.2; 0.3 0.4]))
    @test suitability(suitlayer, (1, 1), 34325) == 0.1
    @test suitability(suitlayer, (2, 2), 7685) == 0.4
end

@testset "suitability sequences are interpolated over timespans" begin

    @testset "sequence cycling" begin
        @test cyclic(-4, 2) == 2
        @test cyclic(-3, 2) == 1
        @test cyclic(-2, 2) == 2
        @test cyclic(-1, 2) == 1
        @test cyclic(0, 2) == 2
        @test cyclic(1, 2) == 1
        @test cyclic(2, 2) == 2
        @test cyclic(3, 2) == 1
        @test cyclic(4, 2) == 2
        @test cyclic(20, 10) == 10
        @test cyclic(21, 10) == 1
        @test cyclic(27, 10) == 7
    end

    # sequence of layers
    global suitseq = SuitabilitySequence(setup.(([0.1 0.2; 0.3 0.4], [0.5 0.6; 0.7 0.8])), 10)

    global ind = ((1, 1),(1, 2),(2, 1),(2, 2))
    @test suitability.((suitseq,), ind, 10) == suitability.((suitseq,), ind, 20)
    @test suitability.((suitseq,), ind, 16) == suitability.((suitseq,), ind, 14)
    @test suitability.((suitseq,), ind, 19) == suitability.((suitseq,), ind, 11)
    @test suitability.((suitseq,), ind, 5)  == suitability.((suitseq,), ind, 45)
    @test suitability.((suitseq,), ind, 15) == suitability.((suitseq,), ind, 55)

    @testset "suitability returns first frame values at 0.5 through the timespan" begin
        @test suitability(suitseq, (1, 1), 5) == 0.1
        @test suitability(suitseq, (1, 2), 5) == 0.2
        @test suitability(suitseq, (2, 1), 5) == 0.3
        @test suitability(suitseq, (2, 2), 5) == 0.4
    end
    @testset "suitability returns second frame values at 1.5 times through the timespan" begin
        @test suitability(suitseq, (1, 1), 15) == 0.5
        @test suitability(suitseq, (1, 2), 15) == 0.6
        @test suitability(suitseq, (2, 1), 15) == 0.7
        @test suitability(suitseq, (2, 2), 15) == 0.8
    end
    @testset "suitability interpolates halfway between frames on the timespan" begin
        @test suitability(suitseq, (1, 1), 10) ≈ 0.3
        @test suitability(suitseq, (1, 2), 10) ≈ 0.4
        @test suitability(suitseq, (2, 1), 10) ≈ 0.5
        @test suitability(suitseq, (2, 2), 10) ≈ 0.6
    end


    global init = [0 0; 0 1]
    global hood = DispersalNeighborhood(; radius=1)
    global model = Models(InwardsBinaryDispersal(neighborhood=hood, prob_threshold=0.0), SuitabilityMask(suitability_threshold=0.4))
    global output = ArrayOutput(init)

    @test Dispersal.pressure(model.models[1], init, 1) 
    Cellular.rule(model.models[1], 0, 2, 2, 2, init, [], suitseq)

    sim!(output, model, init, (suitseq,); time=25)

    # All offset by one, because 1 = t0
    output
    @test output[1]  == setup([0 0; 0 1])
    @test output[2]  == setup([0 0; 1 1])
    @test output[5]  == setup([0 0; 0 1])
    @test output[8]  == setup([0 0; 1 1])
    @test output[10] == setup([0 1; 1 1])
    @test output[15] == setup([1 1; 1 1])
    @test output[20] == setup([0 1; 1 1])
    @test output[22] == setup([0 0; 1 1])
    @test output[25] == setup([0 0; 0 1])
    @test_throws BoundsError output[26]

end

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
        global cc = neighbors(hood, nothing, state, 3, 3, t, source, [])
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

@testset "human dispersal models work" begin

    @testset "Jump dispersal spread randomly" begin
        global layers = SuitabilityLayer(suit)
        global model = Models(JumpDispersal(prob_threshold=0.0, spotrange=3))
        global output = ArrayOutput(init)
        sim!(output, model, init, layers; time=20)
    end

end

@testset "growth models" begin

    global init =  setup([1.0 4.0 7.0;
                          2.0 5.0 8.0;
                          3.0 6.0 9.0])
    
    global output = ArrayOutput(init)


    @testset "fixed growth" begin
        global model = Models(FixedGrowth(2.0))
        sim!(output, model, init; time=3)
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

    @testset "suitability growth" begin
        
        global init =  setup([1.0 1.0 1.0;
                              1.0 1.0 1.0;
                              1.0 1.0 1.0])
        global suit =  setup([1.0 4.0 7.0;
                              2.0 5.0 8.0;
                              3.0 6.0 9.0])
        global layers = SuitabilityLayer(suit)
        global model = Models(SuitabilityGrowth(0.0, 40.0, 1.0))

        sim!(output, model, init, layers; time=3)
        @test output[1] == [1.0 1.0 1.0;
                            1.0 1.0 1.0;
                            1.0 1.0 1.0]
        @test output[2] ≈  [1.0 4.0 7.0;
                            2.0 5.0 8.0;
                            3.0 6.0 9.0]
        @test output[3] ≈  [1.0 16.0 40.0;
                            4.0 25.0 40.0;
                            9.0 36.0 40.0]
    end

end

@testset "jump dispersal models work" begin

    global init =  setup([0 0 0 0 0;
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

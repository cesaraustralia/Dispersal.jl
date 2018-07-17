using Revise, Dispersal, Cellular
using Dispersal: suitability, cyclic, sequence_interpolate, neighbors

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end


@testset "suitability is 1.0 by default" begin
    @test suitability(nothing, (1, 1), 10) == 1.0
    @test suitability(HumanLayer([1 2]), (1, 1), 10) == 1.0
end

@testset "single layer suitability just returns the layer value" begin
    suitlayer = SuitabilityLayer([0.1 0.2; 0.3 0.4])
    @test suitability(suitlayer, (1, 1), 34325) == 0.1
    @test suitability(suitlayer, (2, 2), 7685) == 0.4
end

@testset "suitability sequences are interpolated over timespans" begin
    # sequence of layers
    seq = [[0.1 0.2; 0.3 0.4], 
           [0.5 0.6; 0.7 0.8]]
    suitseq = SuitabilitySequence(seq, 10)

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

    @test length(suitseq.data) == 2
    indices = ((1, 1),(1, 2),(2, 1),(2, 2))
    @test suitability.((suitseq,), indices, 10) == suitability.((suitseq,), indices, 20)
    @test suitability.((suitseq,), indices, 16) == suitability.((suitseq,), indices, 14)
    @test suitability.((suitseq,), indices, 19) == suitability.((suitseq,), indices, 11)
    @test suitability.((suitseq,), indices, 5) == suitability.((suitseq,), indices, 45)
    @test suitability.((suitseq,), indices, 15) == suitability.((suitseq,), indices, 55)

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

    init = [0 0; 0 1] 
    hood = DispersalNeighborhood(; radius=1)
    model = InwardsLocalDispersal(layers=suitseq, neighborhood=hood, prob_threshold=0.0, suitability_threshold=0.4)
    output = ArrayOutput(init)
    sim!(output, model, init; time = 1:30)

    # All offset by one, because 1 = t0
    @test output[1]  == [0 0; 0 1]  
    @test output[2]  == [0 0; 1 1]  
    @test output[6]  == [0 0; 0 1]  
    @test output[9]  == [0 0; 1 1]  
    @test output[11] == [0 1; 1 1]  
    @test output[16] == [1 1; 1 1]  
    @test output[21] == [0 1; 1 1]  
    @test output[23] == [0 0; 1 1]  
    @test output[26] == [0 0; 0 1]  
end

@testset "dispersal kernel array matches passed in function" begin
    f(d) = e^-d
    dk = DispersalNeighborhood(f=f, cellsize=1, radius=2).kernel
    @test typeof(dk) == Array{Float64,2}
    @test size(dk) == (5, 5)

    tesst = [0.0 0.0 0.0 0.0 0.0;
             0.0 0.0 0.0 0.0 0.0;
             0.0 0.0 1.0 0.0 0.0;
             0.0 0.0 0.0 0.0 0.0;
             0.0 0.0 0.0 0.0 0.0;]

    @test sum(dk) ≈ 1.0
    @test_broken dk ≈ test
end

suit =  [1 0 1 1 0;
         0 0 1 1 1;
         1 1 1 1 0;
         1 1 0 1 1;
         1 0 1 1 1]

@testset "dispersal nieghborhood sum matches the passed-in kernel function" begin
    f(d) = 1/d
    hood = DispersalNeighborhood(f=f, cellsize=1, radius=2)
    state = 0
    t = 0

    source = [1 0 0 0 0;
              0 0 0 1 0;
              0 1 0 0 1;
              0 0 0 1 0;
              0 1 0 0 0]

    @testset "neighbouhood sum matches grid*kernel sum for same-sized grid" begin
        cc = neighbors(hood, state, (3, 3), t, source, [])
        @test cc ≈ sum(source .* hood.kernel) 
    end
end

@testset "simple local dispersal simulation" begin

    init =  [0 0 0 0 0;
             0 0 0 0 0;
             0 0 1 0 0;
             0 0 0 0 0;
             0 0 0 0 0]

    test1 = [0 0 0 0 0;
             0 0 0 0 0;
             0 0 1 0 0;
             0 0 0 0 0;
             0 0 0 0 0]

    test2 = [0 0 0 0 0;
             0 0 1 1 0;
             0 1 1 1 0;
             0 1 0 1 0;
             0 0 0 0 0]

    test3 = [0 0 1 1 0;
             0 0 1 1 1;
             1 1 1 1 0;
             1 1 0 1 1;
             1 0 1 1 1]

    # Dispersal in radius 1 neighborhood
    hood = DispersalNeighborhood(; radius=1)
    layers = SuitabilityLayer(suit)

    @testset "inwards dispersal fills the grid where reachable and suitable" begin
        model = InwardsLocalDispersal(layers=layers, neighborhood=hood, prob_threshold=0.0)
        output = ArrayOutput(init)
        sim!(output, model, init; time = 1:3)
        @test output[1] == test1
        @test output[2] == test2
        @test output[3] == test3
    end

    @testset "outwards dispersal fills the grid where reachable and suitable" begin
        model = OutwardsLocalDispersal(layers=layers, neighborhood=hood, prob_threshold=0.0)
        output = ArrayOutput(init)
        sim!(output, model, init; time = 1:3)
        @test output[1] == test1
        @test output[2] == test2
        @test output[3] == test3
    end
end


@testset "jump dispersal models work" begin
    init =  [0 0 0 0 0;
             0 0 0 0 0;
             0 0 1 0 0;
             0 0 0 0 0;
             0 0 0 0 0]

    @testset "Jump dispersal spread randomly" begin
        layers = SuitabilityLayer(suit)
        srand(1234)
        model = JumpDispersal(layers=layers, prob_threshold=0.5)
        output = ArrayOutput(init)
        sim!(output, model, init; time = 1:3)
    end

    @testset "Human dispersal relies on source ans sink population" begin
        layers = SuitabilityLayer(suit)
        srand(1234)
        model = HumanDispersal(layers=layers, prob_threshold=0.5)
        output = ArrayOutput(init)
        sim!(output, model, init; time = 1:3)
    end
end

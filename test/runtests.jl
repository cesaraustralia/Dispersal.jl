using Revise, Dispersal, Cellular
using Dispersal: suitability, cyclic, num_frames

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

@testset "suitability layer interpolation" begin
    # default is 1.0
    @test suitability(nothing, (1, 1), 30) == 1.0
    # sequence of layers
    suitseq = SuitabilitySequence(30, [[0.1 0.2; 0.3 0.4], [0.5 0.6; 0.7 0.8]])
    @test num_frames(suitseq.data) == 2
    @test suitability(suitseq, (1, 1), 15) == 0.1
    @test suitability(suitseq, (2, 2), 15) == 0.4
    @test suitability(suitseq, (1, 2), 45) == 0.6
    @test suitability(suitseq, (2, 1), 45) == 0.7
    # interpolated frames
    @test suitability(suitseq, (2, 2), 30) ≈ 0.6
    @test suitability(suitseq, (1, 1), 30) ≈ 0.3
    suitlayer = SuitabilityLayer([0.1 0.2; 0.3 0.4])
    # single layer suitability
    @test suitability(suitlayer, (1, 1), 34325) == 0.1
    @test suitability(suitlayer, (2, 2), 7685) == 0.4
end

@testset "sequence cycling" begin
    @test cyclic(13, 12) == 1
    @test cyclic(0, 12) == 12
    @test cyclic(3, 12) == 3
end

@testset "build dispersal kernel" begin
    dk = DispersalNeighborhood(f=d->e^-d, radius=1).kernel
    @test typeof(dk) == Array{Float64,2}
    @test size(dk, 1) == 3
    @test size(dk, 2) == 3
    @test dk[1,1] == dk[3,3] == dk[3,1] == dk[1,3]
    @test dk[2,1] == dk[1,2] == dk[3,2] == dk[2,3]
end

@testset "simple local dispersal simulation" begin

    suit =  [1 0 1 1 0;
             0 0 1 1 1;
             1 1 1 1 0;
             1 1 0 1 1;
             1 0 1 1 1]

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

    @testset "inwards" begin
        model = InwardsLocalDispersal(layers=layers, neighborhood=hood, prob=0.0)
        output = ArrayOutput(init)
        sim!(output, model, init; time = 1:3)
        @test output.frames[1] == test1
        @test output.frames[2] == test2
        @test output.frames[3] == test3
    end

    @testset "outwards" begin
        model = OutwardsLocalDispersal(layers=layers, neighborhood=hood, prob=0.0)
        output = ArrayOutput(init)
        sim!(output, model, init; time = 1:3)
        @test output.frames[1] == test1
        @test output.frames[2] == test2
        @test output.frames[3] == test3
    end
end

using Cellular, Dispersal, Test
using Dispersal: get_layers, cyclic

@testset "layers returns 1.0 for an empty tuple" begin
    @test get_layers(nothing, (), (1, 1), 10) == 1
end

suitlayer1 = [0.1 0.2; 
              0.3 0.4]
suitlayer2 = [0.5 0.6; 
              0.7 0.8]

@testset "single layer just returns the layer value" begin
    @test get_layers(nothing, suitlayer1, (1, 1), 34325) == 0.1
    @test get_layers(nothing, suitlayer1, (2, 2), 7685) == 0.4
end

@testset "multiple layers are combined (product)" begin
    @test get_layers(nothing, (suitlayer1, suitlayer2), (1, 1), 1) == 0.1 * 0.5
    @test get_layers(nothing, (suitlayer2, suitlayer2, suitlayer2), (2, 2), 1) == 0.8^3
end


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

suitseq = Sequence(([0.1 0.2; 0.3 0.4], [0.5 0.6; 0.7 0.8]), 10);

@testset "test basic AbstractArray extensions" begin
    @test firstindex(suitseq) == 1
    @test lastindex(suitseq) == 2
    @test size(suitseq) == 2
    @test length(suitseq) == 2
end

@testset "sequence of layers" begin
    timestep = 1
    data = Cellular.FrameData(nothing, nothing, nothing, nothing, timestep, nothing)

    @testset "layers sequences are interpolated over timespans" begin
        ind = ((1, 1),(1, 2),(2, 1),(2, 2))
        @test all(get_layers.(Ref(data), Ref(suitseq), ind, 10) .≈ get_layers.(Ref(data), Ref(suitseq), ind, 20))
        @test all(get_layers.(Ref(data), Ref(suitseq), ind, 16) .≈ get_layers.(Ref(data), Ref(suitseq), ind, 14))
        @test all(get_layers.(Ref(data), Ref(suitseq), ind, 19) .≈ get_layers.(Ref(data), Ref(suitseq), ind, 11))
        @test all(get_layers.(Ref(data), Ref(suitseq), ind, 5)  .≈ get_layers.(Ref(data), Ref(suitseq), ind, 45))
        @test all(get_layers.(Ref(data), Ref(suitseq), ind, 15) .≈ get_layers.(Ref(data), Ref(suitseq), ind, 55))
    end

    @testset "layers returns first frame values at 0.5 through the timespan" begin
        @test get_layers(data, suitseq, (1, 1), 5) == 0.1
        @test get_layers(data, suitseq, (1, 2), 5) == 0.2
        @test get_layers(data, suitseq, (2, 1), 5) == 0.3
        @test get_layers(data, suitseq, (2, 2), 5) == 0.4
    end

    @testset "layers returns second frame values at 1.5 times through the timespan" begin
        @test get_layers(data, suitseq, (1, 1), 15) == 0.5
        @test get_layers(data, suitseq, (1, 2), 15) == 0.6
        @test get_layers(data, suitseq, (2, 1), 15) == 0.7
        @test get_layers(data, suitseq, (2, 2), 15) == 0.8
    end

    @testset "layers interpolates halfway between frames on the timespan" begin
        @test get_layers(data, suitseq, (1, 1), 10) ≈ 0.3
        @test get_layers(data, suitseq, (1, 2), 10) ≈ 0.4
        @test get_layers(data, suitseq, (2, 1), 10) ≈ 0.5
        @test get_layers(data, suitseq, (2, 2), 10) ≈ 0.6
    end
end

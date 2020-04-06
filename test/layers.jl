using Dispersal, Test, Unitful, DimensionalData
using Dispersal: layer, cyclic_index, precalc_timeindex, timestep
using Unitful: d

suitlayer1 = [0.1 0.2;
              0.3 0.4]

@testset "layer returns the layer value" begin
    @test layer(suitlayer1, nothing, (1, 1), 34325) == 0.1
    @test layer(suitlayer1, nothing, (2, 2), 7685) == 0.4
end

@testset "sequence cycling" begin
    @test cyclic_index(-4, 2) == 2
    @test cyclic_index(-3, 2) == 1
    @test cyclic_index(-2, 2) == 2
    @test cyclic_index(-1, 2) == 1
    @test cyclic_index(0, 2) == 2
    @test cyclic_index(1, 2) == 1
    @test cyclic_index(2, 2) == 2
    @test cyclic_index(3, 2) == 1
    @test cyclic_index(4, 2) == 2
    @test cyclic_index(20, 10) == 10
    @test cyclic_index(21, 10) == 1
    @test cyclic_index(27, 10) == 7
end

a = cat([0.1 0.2; 0.3 0.4], [1.5 1.6; 1.7 1.8]; dims=3)
dimz = X(1:2), Y(1:2), Ti(0d:10d:10d)
suitseq = DimensionalArray(a, dimz)

@testset "sequence of layers" begin
    rule = ExactExponentialGrowth(; intrinsicrate=1.0)
    ruleset = Ruleset(rule; timestep=1d)
    data = DynamicGrids.SimData([], ruleset, 0d)

    @testset "layers sequences are interpolated over timespans" begin
        ind = ((1, 1), (1, 2), (2, 1), (2, 2))
        times = (0d, 20d), (16d, 14d), (19d, 11d), (5d, 45d)#, (15d, 55d)

        for (t1, t2) in times
            interp1 = precalc_timeindex(suitseq, rule, data, t1)
            interp2 = precalc_timeindex(suitseq, rule, data, t2)
            @test interp1 == interp2
        end

    end

    @testset "layers returns first frame values at 0.5 through the timespan" begin
        t = 5d
        interp = precalc_timeindex(suitseq, rule, data, t)
        @test layer(suitseq, data, (1, 1), interp) == 0.1
        @test layer(suitseq, data, (1, 2), interp) == 0.2
        @test layer(suitseq, data, (2, 1), interp) == 0.3
        @test layer(suitseq, data, (2, 2), interp) == 0.4
    end

    @testset "layers returns second frame values at 1.5 times through the timespan" begin
        t = 15d
        interp = precalc_timeindex(suitseq, rule, data, t)
        @test layer(suitseq, data, (1, 1), interp) == 1.5
        @test layer(suitseq, data, (1, 2), interp) == 1.6
        @test layer(suitseq, data, (2, 1), interp) == 1.7
        @test layer(suitseq, data, (2, 2), interp) == 1.8
    end

    @testset "layers returns the second frame one timestep in" begin
        t = 10d
        interp = precalc_timeindex(suitseq, rule, data, t)
        @test layer(suitseq, data, (1, 1), interp) ≈ 1.5
        @test layer(suitseq, data, (1, 2), interp) ≈ 1.6
        @test layer(suitseq, data, (2, 1), interp) ≈ 1.7
        @test layer(suitseq, data, (2, 2), interp) ≈ 1.8
    end
end

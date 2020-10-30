using Dispersal, Test, Unitful, DimensionalData
using Dispersal: layer, cyclic_index, precalc_timeindex, timestep
using DynamicGrids: Extent
using Unitful: d

suitlayer1 = [0.1 0.2;
              0.3 0.4]

@testset "layer returns the layer value" begin
    @test layer(suitlayer1, (1, 1), 34325) == 0.1
    @test layer(suitlayer1, (2, 2), 7685) == 0.4
end

@testset "sequence cycling" begin
    @test cyclic_index(-4, 2) == 2
    @test cyclic_index(-3, 2) == 1
    @test cyclic_index(-2, 2) == 2
    @test cyclic_index(-1, 2) == 1
    @test cyclic_index(0, 2)  == 2
    @test cyclic_index(1, 2)  == 1
    @test cyclic_index(2, 2)  == 2
    @test cyclic_index(3, 2)  == 1
    @test cyclic_index(4, 2)  == 2
    @test cyclic_index(-4, 3) == 2
    @test cyclic_index(-3, 3) == 3
    @test cyclic_index(-2, 3) == 1
    @test cyclic_index(-1, 3) == 2
    @test cyclic_index(0, 3)  == 3
    @test cyclic_index(1, 3)  == 1
    @test cyclic_index(2, 3)  == 2
    @test cyclic_index(3, 3)  == 3
    @test cyclic_index(4, 3)  == 1
    @test cyclic_index(20, 10) == 10
    @test cyclic_index(21, 10) == 1
    @test cyclic_index(27, 10) == 7
end



a2 = cat([0.1 0.2; 0.3 0.4], [1.5 1.6; 1.7 1.8]; dims=3)
dimz2 = X(1:2), Y(1:2), Ti(0d:10d:10d)
suitseq2 = DimensionalArray(a2, dimz2)
init = zero(suitseq2[Ti(1)])

a3 = cat([0.1 0.2; 0.3 0.4], [1.5 1.6; 1.7 1.8], [2.1 2.2; 2.3 2.4]; dims=3)
dimz3 = X(1:2), Y(1:2), Ti(0d:10d:20d)
suitseq3 = DimensionalArray(a3, dimz3)
init = zero(suitseq3[Ti(1)])

@testset "layers sequences are interpolated over timespans" begin
    rule = ExactExponentialGrowth(; intrinsicrate=1.0)
    times = (0d, -20d), (16d, 14d), (19d, 11d), (5d, -45d)#, (15d, 55d)
    for (t1, t2) in times
        interp1 = precalc_timeindex(suitseq2, rule, data, t1)
        interp2 = precalc_timeindex(suitseq2, rule, data, t2)
        @test interp1 == interp2
    end

    times = (10d, -20d, 2), (20d, -10d, 3), (0d, 30d, 1) 
    for (t1, t2, v) in times
        interp1 = precalc_timeindex(suitseq3, rule, data, t1)
        interp2 = precalc_timeindex(suitseq3, rule, data, t2)
        @test interp1 == interp2 == v
    end

end

@testset "sequence of layers" begin
    rule = ExactExponentialGrowth(; intrinsicrate=1.0)
    ruleset = Ruleset(rule; timestep=1d)
    data = DynamicGrids.SimData(Extent(init=init, aux=(suit=suitseq3,), tspan=1:1), ruleset)

    @testset "layers returns first frame values at 0.5 through the timespan" begin
        t = 5d
        interp = precalc_timeindex(suitseq3, rule, data, t)
        @test layer(suitseq3, (1, 1), interp) == 0.1
        @test layer(suitseq3, (1, 2), interp) == 0.2
        @test layer(suitseq3, (2, 1), interp) == 0.3
        @test layer(suitseq3, (2, 2), interp) == 0.4
    end

    @testset "layers returns second frame values at 1.5 times through the timespan" begin
        t = 15d
        interp = precalc_timeindex(suitseq3, rule, data, t)
        @test layer(suitseq3, (1, 1), interp) == 1.5
        @test layer(suitseq3, (1, 2), interp) == 1.6
        @test layer(suitseq3, (2, 1), interp) == 1.7
        @test layer(suitseq3, (2, 2), interp) == 1.8
    end

    @testset "layers returns the second frame one timestep in" begin
        t = 10d
        interp = precalc_timeindex(suitseq3, rule, data, t)
        @test layer(suitseq3, (1, 1), interp) ≈ 1.5
        @test layer(suitseq3, (1, 2), interp) ≈ 1.6
        @test layer(suitseq3, (2, 1), interp) ≈ 1.7
        @test layer(suitseq3, (2, 2), interp) ≈ 1.8
    end

end

@testset "start and stop times" begin
    A = rand(4, 5, 10)
    @test Dispersal.timestep(A) == 1
    @test Dispersal.starttime(A) == 1
    @test Dispersal.stoptime(A) == 10
    A = DimensionalArray(rand(4, 5, 10), (X, Y, Ti(20:10:110)))
    @test Dispersal.timestep(A) == 10
    @test Dispersal.starttime(A) == 20
    @test Dispersal.stoptime(A) == 110
end

@testset "LayerCopy" begin
    init = [0 0]
    l = cat([1 2], [3 4]; dims=3)
    @test LayerCopy(:l, 1) === LayerCopy(layerkey=:l)
    @test LayerCopy{:a,:b}(layerkey=:l, timeindex=2) === LayerCopy{:a,:b}(:l, 2)
    ruleset = Ruleset(LayerCopy(layerkey=:l))
    output = ArrayOutput(init; tspan=1:3, aux=(l=l,))
    sim!(output, ruleset)
    @test output == [[0 0], [3 4], [1 2]]
end

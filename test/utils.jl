using Dispersal, Test, Unitful, DimensionalData
using Dispersal: auxval, cyclic_index, precalc_auxtimeindex, timestep
using DynamicGrids: Extent, SimData
using Unitful: d

suitlayer1 = [0.1 0.2;
              0.3 0.4]

@testset "layer returns the layer value" begin
    @test auxval(suitlayer1, 1, 1, 34325) == 0.1
    @test auxval(suitlayer1, 2, 2, 7685) == 0.4
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
init = zero(suitseq)

@testset "sequence of layers" begin
    rule = ExponentialGrowthMap(; auxkey=Val(:suit))
    ruleset = Ruleset(rule; timestep=1d)
    data = SimData(Extent(init=init, aux=(suit=suitseq,), tspan=1:1), ruleset)

    @testset "layers sequences are interpolated over timespans" begin
        ind = ((1, 1), (1, 2), (2, 1), (2, 2))
        times = (0d, 20d), (16d, 14d), (19d, 11d), (5d, 45d)#, (15d, 55d)
        for (t1, t2) in times
            interp1 = precalc_auxtimeindex(suitseq, rule, data, t1)
            interp2 = precalc_auxtimeindex(suitseq, rule, data, t2)
            @test interp1 == interp2
        end
    end

    @testset "auxval returns first frame values at 0.5 through the timespan" begin
        t = 5d
        rule1 = precalc_auxtimeindex(suitseq, rule, data, t)
        t = rule1.auxtimeindex
        @test auxval(suitseq, 1, 1, t) == 0.1
        @test auxval(suitseq, 1, 2, t) == 0.2
        @test auxval(suitseq, 2, 1, t) == 0.3
        @test auxval(suitseq, 2, 2, t) == 0.4
    end

    @testset "auxval returns second frame values at 1.5 times through the timespan" begin
        t = 15d
        rule1 = precalc_auxtimeindex(suitseq, rule, data, t)
        t = rule1.auxtimeindex
        @test auxval(suitseq, 1, 1, t) == 1.5
        @test auxval(suitseq, 1, 2, t) == 1.6
        @test auxval(suitseq, 2, 1, t) == 1.7
        @test auxval(suitseq, 2, 2, t) == 1.8
    end

    @testset "auxval returns the second frame one timestep in" begin
        t = 10d
        rule1 = precalc_auxtimeindex(suitseq, rule, data, t)
        t = rule1.auxtimeindex
        @test auxval(suitseq, 1, 1, t) ≈ 1.5
        @test auxval(suitseq, 1, 2, t) ≈ 1.6
        @test auxval(suitseq, 2, 1, t) ≈ 1.7
        @test auxval(suitseq, 2, 2, t) ≈ 1.8
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

@testset "AuxCopy" begin
    init = [0 0]
    l = cat([1 2], [3 4]; dims=3)
    @test AuxCopy(:l, 1) === AuxCopy(auxkey=:l)
    @test AuxCopy{:a,:b}(auxkey=:l, auxtimeindex=2) === AuxCopy{:a,:b}(:l, 2)
    ruleset = Ruleset(AuxCopy(auxkey=:l))
    output = ArrayOutput(init; tspan=1:3, aux=(l=l,))
    sim!(output, ruleset)
    @test output == [[0 0], [3 4], [1 2]]
end


@testset "Layer on intrinsic growth rate" begin
    
    popSizeInit = [ 1.0 4.0 7.0;
                    2.0 5.0 8.0;
                    3.0 6.0 9.0]

    intrinsicRate = cat([ 1.0 1.0 1.0;
                        1.0 1.0 1.0;
                        1.0 1.0 1.0],
                        [ 2.0 2.0 2.0;
                        2.0 2.0 2.0;
                        2.0 2.0 2.0],
                        [ 1.0 1.0 1.0;
                        1.0 1.0 1.0;
                        1.0 1.0 1.0]; dims=3)

    popSizeGrids = ArrayOutput(popSizeInit; tspan=1:6, aux=(intrinsicRate=intrinsicRate,));
    growthRule = Ruleset(DiscreteGrowthMap(layerkey=:intrinsicRate));
    sim!(popSizeGrids, growthRule);
    @test popSizeGrids[1] == intrinsicRate[:,:,1] .* popSizeInit
    @test popSizeGrids[2] == intrinsicRate[:,:,2] .* popSizeGrids[1]
    @test popSizeGrids[3] == intrinsicRate[:,:,3] .* popSizeGrids[2]
    @test popSizeGrids[4] == intrinsicRate[:,:,1] .* popSizeGrids[3]
    @test popSizeGrids[5] == intrinsicRate[:,:,2] .* popSizeGrids[4]
    @test popSizeGrids[6] == intrinsicRate[:,:,3] .* popSizeGrids[5]


    popSizeInit = [ 1.0 4.0 7.0;
                    2.0 5.0 8.0;
                    3.0 6.0 9.0]

    intrinsicRate = cat([ 1.0 1.0 1.0;
                        1.0 1.0 1.0;
                        1.0 1.0 1.0],
                        [ 2.0 2.0 2.0;
                        2.0 2.0 2.0;
                        2.0 2.0 2.0],
                        [ 1.0 1.0 1.0;
                        1.0 1.0 1.0;
                        1.0 1.0 1.0]; dims=3)

    popSizeGrids = ArrayOutput(popSizeInit; tspan=1:6, aux=(intrinsicRate=intrinsicRate,));
    growthRule = Ruleset(ExactExponentialGrowthMap(layerkey=:intrinsicRate));
    sim!(popSizeGrids, growthRule);
    @test popSizeGrids[1] == popSizeInit
    @test popSizeGrids[2] == popSizeGrids[1].*exp.(intrinsicRate[:,:,2])
    @test popSizeGrids[3] == popSizeGrids[2].*exp.(intrinsicRate[:,:,3])
    @test popSizeGrids[4] == popSizeGrids[3].*exp.(intrinsicRate[:,:,1])
    @test popSizeGrids[5] == popSizeGrids[4].*exp.(intrinsicRate[:,:,2])
    @test popSizeGrids[6] == popSizeGrids[5].*exp.(intrinsicRate[:,:,3])

end

@testset "Test double layer on intrinsic growth rate and carying capacity" begin
    popSizeInit = [ 1.0 4.0 7.0;
                    2.0 5.0 8.0;
                    3.0 6.0 9.0]

    intrinsicRate = cat([ 1.0 1.0 1.0;
                        1.0 1.0 1.0;
                        1.0 1.0 1.0],
                        [ 2.0 2.0 2.0;
                        2.0 2.0 2.0;
                        2.0 2.0 2.0],
                        [ 1.0 1.0 1.0;
                        1.0 1.0 1.0;
                        1.0 1.0 1.0]; dims=3)

    carryingCapacity = cat([ 10.0 10.0 10.0;
                            10.0 10.0 10.0;
                            10.0 10.0 10.0],
                            [ 10.0 10.0 10.0;
                            10.0 10.0 10.0;
                            10.0 10.0 10.0],
                            [ 10.0 10.0 10.0;
                            10.0 10.0 10.0;
                            10.0 10.0 10.0]; dims=3)

    popParameter = cat(intrinsicRate, carryingCapacity; dims = 4)

    popSizeGrids = ArrayOutput(popSizeInit; tspan=1:6, aux=(popParameter=popParameter,));
    growthRule = Ruleset(ExactLogisticGrowthMap2(layerkey=:popParameter));
    sim!(popSizeGrids, growthRule);
end
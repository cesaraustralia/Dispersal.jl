using Dispersal

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# write your own tests here
@testset "suitability layer interpolation" begin
    # default is 1.0
    @test Dispersal.suitability(nothing, 1, 1, 30) == 1.0
    # sequence of layers
    suitseq = SuitabilitySequence(30, [[0.1 0.2; 0.3 0.4], [0.5 0.6; 0.7 0.8]])
    @test Dispersal.frames(suitseq.data) == 2
    @test Dispersal.suitability(suitseq, 1, 1, 15) == 0.1
    @test Dispersal.suitability(suitseq, 2, 2, 15) == 0.4
    @test Dispersal.suitability(suitseq, 1, 2, 45) == 0.6
    @test Dispersal.suitability(suitseq, 2, 1, 45) == 0.7
    # interpolated frames
    @test Dispersal.suitability(suitseq, 2, 2, 30) ≈ 0.6
    @test Dispersal.suitability(suitseq, 1, 1, 30) ≈ 0.3
    suitlayer = SuitabilityLayer([0.1 0.2; 0.3 0.4])
    # single layer suitability
    @test Dispersal.suitability(suitlayer, 1, 1, 34325) == 0.1
    @test Dispersal.suitability(suitlayer, 2, 2, 7685) == 0.4
end

@test Dispersal.cyclic(13, 12) == 1
@test Dispersal.cyclic(0, 12) == 12

@testset "build dispersal kernel" begin
    dk = Dispersal.build_dispersal_kernel(d->e^-d, 1)
    @test typeof(dk) == Array{Float64,2}
    @test size(dk, 1) == 3
    @test size(dk, 2) == 3
    @test dk[1,1] == dk[3,3] == dk[3,1] == dk[1,3]
    @test dk[2,1] == dk[1,2] == dk[3,2] == dk[2,3]
end


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
    Cellular.rule(model.models[1], 0, (2, 2), 2, init, [], suitseq)

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

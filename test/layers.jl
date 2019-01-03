using Dispersal: get_layers, cyclic

@testset "layers returns 1.0 by default" begin
    @test get_layers(nothing, (1, 1), 10) == 1.0
end

@testset "single layer just returns the layer value" begin
    global suitlayer = setup([0.1 0.2; 0.3 0.4])
    @test get_layers(suitlayer, (1, 1), 34325) == 0.1
    @test get_layers(suitlayer, (2, 2), 7685) == 0.4
end

@testset "layers sequences are interpolated over timespans" begin

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
    global suitseq = Sequence(setup.(([0.1 0.2; 
                                       0.3 0.4], 
                                      [0.5 0.6; 
                                       0.7 0.8])), 
                              10);

    global ind = ((1, 1),(1, 2),(2, 1),(2, 2))
    @test get_layers.((suitseq,), ind, 10) == get_layers.((suitseq,), ind, 20)
    @test get_layers.((suitseq,), ind, 16) == get_layers.((suitseq,), ind, 14)
    @test get_layers.((suitseq,), ind, 19) == get_layers.((suitseq,), ind, 11)
    @test get_layers.((suitseq,), ind, 5)  == get_layers.((suitseq,), ind, 45)
    @test get_layers.((suitseq,), ind, 15) == get_layers.((suitseq,), ind, 55)

    @testset "layers returns first frame values at 0.5 through the timespan" begin
        @test get_layers(suitseq, (1, 1), 5) == 0.1
        @test get_layers(suitseq, (1, 2), 5) == 0.2
        @test get_layers(suitseq, (2, 1), 5) == 0.3
        @test get_layers(suitseq, (2, 2), 5) == 0.4
    end
    @testset "layers returns second frame values at 1.5 times through the timespan" begin
        @test get_layers(suitseq, (1, 1), 15) == 0.5
        @test get_layers(suitseq, (1, 2), 15) == 0.6
        @test get_layers(suitseq, (2, 1), 15) == 0.7
        @test get_layers(suitseq, (2, 2), 15) == 0.8
    end
    @testset "layers interpolates halfway between frames on the timespan" begin
        @test get_layers(suitseq, (1, 1), 10) ≈ 0.3
        @test get_layers(suitseq, (1, 2), 10) ≈ 0.4
        @test get_layers(suitseq, (2, 1), 10) ≈ 0.5
        @test get_layers(suitseq, (2, 2), 10) ≈ 0.6
    end


    global init = [0 0; 0 1]
    global hood = DispersalKernel(; radius=1)
    global model = Models(InwardsBinaryDispersal(neighborhood=hood, prob_threshold=0.0), 
                          SuitabilityMask(layers=suitseq, threshold=0.4))
    global output = ArrayOutput(init, 25)

    @test Dispersal.pressure(model.models[1], init, 1) 

    sim!(output, model, init; tstop=25)

    # All offset by one, because 1 = t0
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

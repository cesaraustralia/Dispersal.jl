using Dispersal, Test, Unitful, Dates, DimensionalData
using Unitful: d, s

struct TestFormulation <: KernelFormulation end
(f::TestFormulation)(d) = 1.0

@testset "dispersal kernel array matches passed in function" begin
    init = [0.0 1.0 0.0;
            1.0 0.0 1.0]
    radius = 2
    dk = DispersalKernel{radius}(formulation=ExponentialKernel(1.0), cellsize=1.0).kernel
    @test size(dk) == (5, 5)
    @test sum(dk) ≈ 1.0
end

@testset "binary dispersal and growth mask" begin
    init = [0.0 0.0; 0.0 1.0]
    radius = 1
    hood = DispersalKernel{radius}(; formulation=ExponentialKernel(2.0))


    # time sequence for auxillary input
    a = cat([0.1 0.2; 0.3 0.4], [0.5 0.6; 0.7 0.8], dims=3)
    dimz = Y(1:2), X(1:2), Ti(1d:10d:11d)
    suitseq = DimArray(a, dimz)

    # Regular
    ruleset1 = Ruleset(InwardsDispersal(neighborhood=hood),
                       ThresholdGrowth(rate=Aux(:suitseq), threshold=0.4); timestep=1d)
    # Chained
    ruleset2 = Ruleset(Chain(InwardsDispersal(neighborhood=hood),
                       ThresholdGrowth(rate=Aux(:suitseq), threshold=0.4)) ; timestep=1d)
    output1 = ArrayOutput(init; tspan=1d:1d:25d)
    output2 = ArrayOutput(init; tspan=1d:1d:25d)

    sim!(output1, ruleset1; aux=(suitseq=suitseq,))
    sim!(output2, ruleset2; aux=(suitseq=suitseq,))

    results = [[0 0; 0 1],
               [0 0; 0 1],
               [0 0; 0 1],
               [0 0; 0 1],
               [0 0; 0 1],
               [1 1; 1 1],
               [1 1; 1 1],
               [0 0; 0 1],
               [0 0; 0 1],]

    @testset "Rounded" begin
        @test map(o -> o .> 0, output1[[1, 2, 5, 8, 10, 15, 20, 22, 25]]) == results
        @test map(o -> o .> 0, output2[[1, 2, 5, 8, 10, 15, 20, 22, 25]]) == results
        output1[[1, 2, 5, 8, 10, 15, 20, 22, 25]]
        output2[[1, 2, 5, 8, 10, 15, 20, 22, 25]]
    end

    # @testset "Interpolated" begin
    #     # All offset by one, because 1 = t0
    #     @test output1[1]  == output2[1]  == [0 0; 0 1]
    #     @test output1[2]  == output2[2]  == [0 0; 1 1]
    #     @test output1[5]  == output2[5]  == [0 0; 0 1]
    #     @test output1[8]  == output2[8]  == [0 0; 1 1]
    #     @test output1[10] == output2[10] == [0 1; 1 1]
    #     @test output1[15] == output2[15] == [1 1; 1 1]
    #     @test output1[20] == output2[20] == [0 1; 1 1]
    #     @test output1[22] == output2[22] == [0 0; 1 1]
    #     @test output1[25] == output2[25] == [0 0; 0 1]
    # end
    
    @test_throws BoundsError output1[26]
    @test_throws BoundsError output2[26]
end

@testset "floating point population dispersal simulation with suitability mask" begin

    suit =  [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
             0.0 0.0 1.0 0.0 1.0 0.0 1.0 1.0 0.0 0.0 0.0;
             0.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0 1.0 0.0 0.0;
             0.0 0.0 1.0 0.0 1.0 1.0 1.0 1.0 0.0 0.0 0.0;
             0.0 0.0 1.0 0.0 1.0 1.0 0.0 1.0 1.0 0.0 0.0;
             0.0 0.0 1.0 1.0 1.0 0.0 1.0 1.0 1.0 0.0 0.0;
             0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
             0.0 0.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 0.0 0.0;
             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;]

    init =  [0.0 0.0   0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0 0.0;
             0.0 0.0   0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0 0.0;
             0.0 0.0 100.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0 0.0;
             0.0 0.0   0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0 0.0;
             0.0 0.0   0.0 0.0 0.0 0.0 100.0 0.0 0.0 0.0 0.0;
             0.0 0.0   0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0 0.0;
             0.0 0.0   0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0 0.0;
             0.0 0.0   0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0 0.0;
             0.0 0.0   0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0 0.0;
             0.0 0.0   0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0 0.0;
             0.0 0.0   0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0 0.0;]

    test1 = [0.0 0.0   0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0 0.0;
             0.0 0.0   0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0 0.0;
             0.0 0.0 100.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0 0.0;
             0.0 0.0   0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0 0.0;
             0.0 0.0   0.0 0.0 0.0 0.0 100.0 0.0 0.0 0.0 0.0;
             0.0 0.0   0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0 0.0;
             0.0 0.0   0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0 0.0;
             0.0 0.0   0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0 0.0;
             0.0 0.0   0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0 0.0;
             0.0 0.0   0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0 0.0;
             0.0 0.0   0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0 0.0;]

    test2 = [0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0;
             0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0;
             0.0  0.0  4.0  0.0  8.0  0.0  4.0  4.0  0.0  0.0  0.0;
             0.0  0.0  4.0  4.0  0.0  0.0  4.0  4.0  4.0  0.0  0.0;
             0.0  0.0  4.0  0.0  8.0  4.0  4.0  4.0  0.0  0.0  0.0;
             0.0  0.0  0.0  0.0  4.0  4.0  0.0  4.0  4.0  0.0  0.0;
             0.0  0.0  0.0  0.0  4.0  0.0  4.0  4.0  4.0  0.0  0.0;
             0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0;
             0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0;
             0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0;
             0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0;]

    test3 = [0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0;
             0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0;
             0.0  0.0  1.28 0.0  1.92 0.0  1.92 1.28 0.0  0.0  0.0;
             0.0  0.0  1.44 1.76 0.0  0.0  2.56 1.76 1.44 0.0  0.0;
             0.0  0.0  1.6  0.0  2.56 2.88 3.2  2.24 0.0  0.0  0.0;
             0.0  0.0  1.12 0.0  1.92 2.24 0.0  1.92 1.6  0.0  0.0;
             0.0  0.0  0.8  1.12 1.44 0.0  2.08 1.44 1.12 0.0  0.0;
             0.0  0.0  0.32 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0;
             0.0  0.0  0.16 0.16 0.32 0.48 0.64 0.48 0.48 0.0  0.0;
             0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0;
             0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0;]

    # Dispersal in radius 1 neighborhood
    mask = ThresholdGrowth(rate=Aux(:suit))
    rad = 2

    @testset "inwards population dispersal fills the grid where reachable and suitable" begin
        hood = DispersalKernel{rad}(; formulation=TestFormulation(), distancemethod=CentroidToCentroid())
        inwards = InwardsDispersal(neighborhood=hood)
        rules = Ruleset(inwards, mask; timestep=1d)
        output = ArrayOutput(init; tspan=4d:1d:6d, aux=(suit=suit,))
        for kw in ((), (opt=SparseOpt(),), (proc=ThreadedCPU(),), (proc=ThreadedCPU(), opt=SparseOpt()))
            sim!(output, rules; kw...)
            @test output[1] == test1
            @test output[2] == test2
            @test output[3] ≈ test3
            # Chained
            rules = Ruleset(Chain(inwards, mask); timestep=1d, kw...)
            sim!(output, rules; tspan=4d:1d:6d)
            @test output[1] == test1
            @test output[2] == test2
            @test output[3] ≈ test3
        end
    end

    @testset "outwards population dispersal fills the grid where reachable and suitable" begin
        hood = DispersalKernel{rad}(; formulation=TestFormulation(), distancemethod=CentroidToCentroid())
        outwards = OutwardsDispersal(neighborhood=hood)
        rules = Ruleset(outwards, mask; timestep=Month(1))
        output = ArrayOutput(init; tspan=Date(2001, 1):Month(1):Date(2001, 3), aux=(suit=suit,))
        for kw in ((), (opt=SparseOpt(),), (proc=ThreadedCPU(),), (proc=ThreadedCPU(), opt=SparseOpt()))
            sim!(output, rules; kw...)
            @test output[1] == test1
            @test output[2] == test2
            @test output[3] ≈ test3
        end
    end

end

@testset "OutwardsDispersal Test" begin
    # Define a mask
    mask_data = [i == 1 || i == 10 || j == 1 || j == 10 ? false : true for i in 1:10, j in 1:10]

    # Create OutwardsDispersal with a mask
    outdisp_with_mask = OutwardsDispersal(
        formulation=ExponentialKernel(λ=0.0125),
        distancemethod=AreaToArea(30),
        mask_flag=Mask()
    )

    # Create OutwardsDispersal without a mask, NoMask is default
    outdisp_without_mask = OutwardsDispersal(
        formulation=ExponentialKernel(λ=0.0125),
        distancemethod=AreaToArea(30)
    )

    # Create a grid with empty borders matching the mask
    init = map(x -> x ? 100.0 : 0.0, mask_data)

    # Create ruleset and outputs
    rule_with_mask = Ruleset(outdisp_with_mask; boundary=Reflect())
    rule_without_mask = Ruleset(outdisp_without_mask; boundary=Reflect())

    # Run the simulation with a mask
    output_with_mask = ArrayOutput(init; tspan=1:1000, mask=mask_data)
    a = sim!(output_with_mask, rule_with_mask)
    
    @test sum(a[1]) ≈ sum(a[1000]) # Floating error should be smaller than 1.0

    # Run the simulation without a mask to check default works fine
    output_without_mask = ArrayOutput(init; tspan=1:1000)
    b = sim!(output_without_mask, rule_without_mask)

    @test sum(b[1]) ≈ sum(b[1000]) # Floating error should be smaller than 1.0

    # Run the simulation with a mask but outdisp_without_mask
    output_without_mask = ArrayOutput(init; tspan=1:1000, mask=mask_data)
    b = sim!(output_without_mask, rule_without_mask)

    @test sum(b[1]) !≈ sum(b[1000]) #= Floating error should be larger than 1.0
    because this does not identify the mask properly =#
end
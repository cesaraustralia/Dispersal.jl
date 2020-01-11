using Dispersal, Test, Unitful, Dates
using DimensionalData: DimensionalArray, X, Y, Time
using Unitful: d, s

struct TestFormulation <: KernelFormulation end
(f::TestFormulation)(d) = 1.0

@testset "dispersal kernel array matches passed in function" begin
    init = [0.0 1.0 0.0;
            1.0 0.0 1.0]
    radius = 2
    dk = DispersalKernel{radius}(formulation=ExponentialKernel(1.0), kernel=init, cellsize=1.0).kernel
    @test size(dk) == (5, 5)
    @test sum(dk) ≈ 1.0
end

@testset "binary dispersal and mask" begin
    init = [0 0; 0 1]
    radius = 1
    hood = DispersalKernel{radius}()

    # time sequence of layers
    a = cat([0.1 0.2; 0.3 0.4], [0.5 0.6; 0.7 0.8], dims=3)
    dimz = X(1:2), Y(1:2), Time(1d:10d:11d)
    suitseq = DimensionalArray(a, dimz)

    # Regular
    ruleset1 = Ruleset(InwardsBinaryDispersal(neighborhood=hood, prob_threshold=0.0),
                       MaskGrowthMap(layer=suitseq, threshold=0.4); init=init, timestep=1d)
    # Chained
    ruleset2 = Ruleset(Chain(InwardsBinaryDispersal(neighborhood=hood, prob_threshold=0.0),
                       MaskGrowthMap(layer=suitseq, threshold=0.4)); init=init, timestep=1d)
    output1 = ArrayOutput(init, 25)
    output2 = ArrayOutput(init, 25)

    sim!(output1, ruleset1; tspan=(1d, 25d))
    sim!(output2, ruleset2; tspan=(1d, 25d))

    @testset "Rounded" begin
        @test output1[1]  == [0 0; 0 1]
        @test output1[2]  == [0 0; 0 1]
        @test output1[5]  == [0 0; 0 1]
        @test output1[8]  == [0 0; 0 1]
        @test output1[10] == [0 0; 0 1]
        @test output1[15] == [1 1; 1 1]
        @test output1[20] == [1 1; 1 1]
        @test output1[22] == [0 0; 0 1]
        @test output1[25] == [0 0; 0 1]
        @test output2[1]  == [0 0; 0 1]
        @test output2[2]  == [0 0; 0 1]
        @test output2[5]  == [0 0; 0 1]
        @test output2[8]  == [0 0; 0 1]
        @test output2[10] == [0 0; 0 1]
        @test output2[15] == [1 1; 1 1]
        @test output2[20] == [1 1; 1 1]
        @test output2[22] == [0 0; 0 1]
        @test output2[25] == [0 0; 0 1]
    end

    # @testset "Interpolated" begin
    #     # All offset by one, because 1 = t0
    #     @test output1[1]  == [0 0; 0 1]
    #     @test output1[2]  == [0 0; 1 1]
    #     @test output1[5]  == [0 0; 0 1]
    #     @test output1[8]  == [0 0; 1 1]
    #     @test output1[10] == [0 1; 1 1]
    #     @test output1[15] == [1 1; 1 1]
    #     @test output1[20] == [0 1; 1 1]
    #     @test output1[22] == [0 0; 1 1]
    #     @test output1[25] == [0 0; 0 1]
    #     # All offset by one, because 1 = t0
    #     @test output2[1]  == [0 0; 0 1]
    #     @test output2[2]  == [0 0; 1 1]
    #     @test output2[5]  == [0 0; 0 1]
    #     @test output2[8]  == [0 0; 1 1]
    #     @test output2[10] == [0 1; 1 1]
    #     @test output2[15] == [1 1; 1 1]
    #     @test output2[20] == [0 1; 1 1]
    #     @test output2[22] == [0 0; 1 1]
    #     @test output2[25] == [0 0; 0 1]
    # end
    
    @test_throws BoundsError output1[26]
    @test_throws BoundsError output2[26]
end


@testset "binary dispersal simulation with suitability mask" begin

    suit =  [1 0 1 1 0;
             0 0 1 1 1;
             1 1 1 1 0;
             1 1 0 1 1;
             1 0 1 1 1]

    init =  Bool[0 0 0 0 0;
                 0 0 0 0 0;
                 0 0 1 0 0;
                 0 0 0 0 0;
                 0 0 0 0 0]

    test1 = Bool[0 0 0 0 0;
                 0 0 0 0 0;
                 0 0 1 0 0;
                 0 0 0 0 0;
                 0 0 0 0 0]

    test2 = Bool[0 0 0 0 0;
                 0 0 1 1 0;
                 0 1 1 1 0;
                 0 1 0 1 0;
                 0 0 0 0 0]

    test3 = Bool[0 0 1 1 0;
                 0 0 1 1 1;
                 1 1 1 1 0;
                 1 1 0 1 1;
                 1 0 1 1 1]

    # Dispersal in radius 1 neighborhood
    mask = MaskGrowthMap(layer=suit)
    radius = 1
    hood = DispersalKernel{radius}(; formulation=ExponentialKernel(1.0))

    @testset "inwards binary dispersal fills the grid where reachable and suitable" begin
        inwards = InwardsBinaryDispersal(neighborhood=hood, prob_threshold=0.0)
        rules = Ruleset(inwards, mask; init=init)
        output = ArrayOutput(init, 3)
        sim!(output, rules; tspan=(1, 3))
        @test output[1] == test1
        @test output[2] == test2
        @test output[3] == test3

        # As subrules
        rules = Ruleset(inwards, mask; init=init)
        sim!(output, rules; tspan = (1, 3))
        @test output[1] == test1
        @test output[2] == test2
        @test output[3] == test3
    end
    @testset "outwards dispersal fills the grid where reachable and mask" begin
        outwards = OutwardsBinaryDispersal(neighborhood=hood, prob_threshold=0.0)
        rules = Ruleset(outwards, mask; init=init)
        output = ArrayOutput(init, 3)
        sim!(output, rules; tspan=(1, 3))
        @test output[1] == test1
        @test output[2] == test2
        @test output[3] == test3
    end

end

@testset "floating point population dispersal simulation with suitability mask" begin

    suit =  [1.0 0.0 1.0 0.0 1.0 1.0 0.0;
             1.0 1.0 0.0 0.0 1.0 1.0 1.0;
             1.0 0.0 1.0 1.0 1.0 1.0 0.0;
             1.0 0.0 1.0 1.0 0.0 1.0 1.0;
             1.0 1.0 1.0 0.0 1.0 1.0 1.0;
             1.0 0.0 0.0 0.0 0.0 0.0 0.0;
             1.0 1.0 1.0 1.0 1.0 1.0 1.0]

    init =  [100.0 0.0 0.0 0.0 0.0 0.0 0.0;
             0.0 0.0 0.0 0.0   0.0 0.0 0.0;
             0.0 0.0 0.0 0.0 100.0 0.0 0.0;
             0.0 0.0 0.0 0.0   0.0 0.0 0.0;
             0.0 0.0 0.0 0.0   0.0 0.0 0.0;
             0.0 0.0 0.0 0.0   0.0 0.0 0.0;
             0.0 0.0 0.0 0.0   0.0 0.0 0.0]

    test1 = [100.0 0.0 0.0 0.0 0.0 0.0 0.0;
             0.0 0.0 0.0 0.0   0.0 0.0 0.0;
             0.0 0.0 0.0 0.0 100.0 0.0 0.0;
             0.0 0.0 0.0 0.0   0.0 0.0 0.0;
             0.0 0.0 0.0 0.0   0.0 0.0 0.0;
             0.0 0.0 0.0 0.0   0.0 0.0 0.0;
             0.0 0.0 0.0 0.0   0.0 0.0 0.0]

    test2 = [4.0  0.0  8.0  0.0  4.0  4.0  0.0;
             4.0  4.0  0.0  0.0  4.0  4.0  4.0;
             4.0  0.0  8.0  4.0  4.0  4.0  0.0;
             0.0  0.0  4.0  4.0  0.0  4.0  4.0;
             0.0  0.0  4.0  0.0  4.0  4.0  4.0;
             0.0  0.0  0.0  0.0  0.0  0.0  0.0;
             0.0  0.0  0.0  0.0  0.0  0.0  0.0]

    test3 = [1.28  0.0   1.92  0.0   1.92  1.28  0.0;
             1.44  1.76  0.0   0.0   2.56  1.76  1.44;
             1.6   0.0   2.56  2.88  3.2   2.24  0.0;
             1.12  0.0   1.92  2.24  0.0   1.92  1.6;
             0.8   1.12  1.44  0.0   2.08  1.44  1.12;
             0.32  0.0   0.0   0.0   0.0   0.0   0.0;
             0.16  0.16  0.32  0.48  0.64  0.48  0.48;]

    # Dispersal in radius 1 neighborhood
    mask = MaskGrowthMap(layer=suit)
    radius = 2

    @testset "inwards population dispersal fills the grid where reachable and suitable" begin
        hood = DispersalKernel{radius}(; formulation=TestFormulation(), distancemethod=CentroidToCentroid())
        inwards = InwardsPopulationDispersal(neighborhood=hood)
        rules = Ruleset(inwards, mask; init=init, timestep=1d)
        output = ArrayOutput(init, 3)
        sim!(output, rules; tspan=(4d, 6d))
        @test output[1] == test1
        @test output[2] == test2
        @test output[3] ≈ test3

        # As subrules
        rules = Ruleset(Chain(inwards, mask); init=init, timestep=Day(1))
        sim!(output, rules; tspan=(Date(2001,1,1), Date(2001,1,3)))
        @test output[1] == test1
        @test output[2] == test2
        @test output[3] ≈ test3
    end

    @testset "outwards population dispersal fills the grid where reachable and suitable" begin
        hood = DispersalKernel{radius}(; formulation=TestFormulation())
        outwards = OutwardsPopulationDispersal(neighborhood=hood)
        rules = Ruleset(outwards, mask; init=init, timestep=Month(1))
        output = ArrayOutput(init, 3)
        sim!(output, rules; tspan=(Date(2001, 1), Date(2001, 3)))
        @test output[1] == test1
        @test output[2] == test2
        @test_broken output[3] ≈ test3
    end
end

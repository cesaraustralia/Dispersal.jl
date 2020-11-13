using Dispersal, Test

@testset "discrete growth" begin
    init =  [1.0 4.0 7.0;
             2.0 5.0 8.0;
             3.0 6.0 9.0]

    test1 = init
    test2 = init .* 2.0
    test3 = init .* 2.0 .^2

    output = ArrayOutput(init; tspan=1:3)
    rule = DiscreteGrowth(intrinsicrate=2.0)
    sim!(output, rule)

    @test output[1] == test1
    @test output[2] == test2
    @test output[3] == test3
end

@testset "discrete growth map" begin
    init =  [1.0 4.0 7.0;
             2.0 5.0 8.0;
             3.0 6.0 9.0]

    suit =  [1.0 1.0 2.0;
             2.0 1.0 0.5;
             1.0 1.0 0.5]

    output = ArrayOutput(init; tspan=1:3, aux=(suit=suit,))
    rule = DiscreteGrowthMap(; ratekey=Val(:suit))
    sim!(output, rule)

    @test output[1] == init
    @test output[2] == init .* suit
    @test output[3] == init .* suit .^2
end

@testset "Test double layers LogisticGrowthMap3" begin
    pop_size_init = [1.0 4.0 7.0;
                     2.0 5.0 8.0;
                     3.0 6.0 9.0]

    intrinsicrate = cat([1.0 1.0 1.0
                         1.0 1.0 1.0
                         1.0 1.0 1.0],
                        [2.0 2.0 2.0
                         2.0 2.0 2.0
                         2.0 2.0 2.0],
                        [1.0 1.0 1.0
                         1.0 1.0 1.0
                         1.0 1.0 1.0]; dims=3)

    carrying_capacity = cat([10.0 10.0 10.0
                             10.0 10.0 10.0
                             10.0 10.0 10.0],
                            [10.0 10.0 10.0
                             10.0 10.0 10.0
                             10.0 10.0 10.0],
                            [10.0 10.0 10.0
                             10.0 10.0 10.0
                             10.0 10.0 10.0]; dims=3)

    output = ArrayOutput(pop_size_init;
        tspan=1:6, aux=(intrinsicrate=intrinsicrate, carrycap=carrying_capacity,)
    )
    rule = Dispersal.LogisticGrowthMap3(;
        ratekey=Val(:intrinsicrate),
        carrycapkey=Val(:carrycap),
    )
    sim!(output, rule);
end

@testset "Test double grids DiscreteGrowth" begin
    init = (
        pop1 = [0.5 0.0 0.0
                0.0 0.0 1.0
                0.0 0.0 0.0],
        pop2 = [0.5 0.0 0.0
                0.0 1.0 0.0
                0.0 0.0 0.0],
    )

    exposure = [1.0 0.5 0.0
                0.5 0.0 0.0
                0.0 0.0 0.0]

    output = ArrayOutput(init; tspan=1:3)

    capped = Dispersal.CappedDiscreteGrowth{Tuple{:pop1,:pop2},Tuple{:pop1,:pop2}}(
        intrinsicrate=(1.5, 2.0), carrycap=10.0
    )
    hood = DispersalKernel{1}(; formulation=ExponentialKernel(15))
    dispersal_pop1 = InwardsPopulationDispersal{:pop1,:pop1}(;neighborhood=hood)
    dispersal_pop2 = InwardsPopulationDispersal{:pop2,:pop2}(;neighborhood=hood)

    ruleset = Ruleset(
       capped, dispersal_pop1, dispersal_pop2;
       overflow=WrapOverflow()
    )

    sim!(output, ruleset)
    output[1].pop1
    output[2].pop1

    output[1].pop2
    output[2].pop2
end

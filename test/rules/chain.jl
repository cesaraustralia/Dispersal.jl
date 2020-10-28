
using Dispersal, Test

@testset "chain growth and allee" begin
    init = [1.0 4.0 7.0;
            2.0 5.0 8.0;
            3.0 6.0 9.0]

    output = ArrayOutput(init; tspan=1:3)

    growthRule = ExactExponentialGrowth(intrinsicrate=-log(2.0))
    alleeRule = AlleeExtinction(minfounders = 2.0)
    ruleset = Ruleset(Chain(growthRule, alleeRule))
    sim!(output, ruleset)
    @test output[1] == init
    @test output[2] == [0.0  2.0  3.5;
                        0.0  2.5  4.0;
                        0.0  3.0  4.5]
    @test output[3] == [0.0  0.0  0.0;
                        0.0  0.0  2.0;
                        0.0  0.0  2.25]
end


@testset "chain discrete growth and survival" begin
    init = [1.0 4.0 7.0;
            2.0 5.0 8.0;
            3.0 6.0 9.0]

    suit = [10.0 10.0 10.0;
            5.0 5.0 5.0;
            2.0 2.0 2.0]

    output = ArrayOutput(init; tspan=1:3,aux=(suit=suit,))

    growthRule = DiscreteGrowth(intrinsicrate=2.0)
    survivalRule = SurvLogLogisticMap(layerkey=Val(:suit), LC50=5.0)
    ruleset = Ruleset(Chain(growthRule, survivalRule))
    sim!(output, ruleset)
    @test output[1] == init
    @test output[2] ≈ [0.965357  3.86143  6.7575;
                        2.0       5.0      8.0;
                        3.13735   6.2747   9.41204] atol=1e-4
    @test output[3] ≈ [0.931913  3.72765  6.52339;
                        2.0       5.0      8.0;
                        3.28098   6.56197  9.84295] atol=1e-4
end

@testset "chain SimData growth and survival" begin
    init = [1.0 4.0 7.0;
            2.0 5.0 8.0;
            3.0 6.0 9.0]

    suit = [10.0 10.0 10.0;
            5.0 5.0 5.0;
            2.0 2.0 2.0]

    output = ArrayOutput(init; tspan=1:3,aux=(suit=suit,))

    growthRule = DiscreteGrowth(intrinsicrate=2.0)
    survivalRule = SurvLogLogisticMap(layerkey=Val(:suit), LC50=5.0)
    ruleset = Ruleset(Chain(growthRule, survivalRule))
    sim!(output, ruleset)
    @test output[1] == init
    @test output[2] ≈ [0.965357  3.86143  6.7575;
                        2.0       5.0      8.0;
                        3.13735   6.2747   9.41204] atol=1e-4
    @test output[3] ≈ [0.931913  3.72765  6.52339;
                        2.0       5.0      8.0;
                        3.28098   6.56197  9.84295] atol=1e-4
end
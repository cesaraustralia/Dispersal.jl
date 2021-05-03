
using Dispersal, Test

@testset "chain growth and allee" begin
    init = [1.0 4.0 7.0;
            2.0 5.0 8.0;
            3.0 6.0 9.0]

    output = ArrayOutput(init; tspan=1:3)

    growth = ExponentialGrowth(; rate=-log(2.0), timestep=1)
    allee = AlleeExtinction(minfounders = 2.0)
    ruleset = Ruleset(Chain(growth, allee))
    sim!(output, ruleset)
    @test output[1] == init
    @test output[2] == [0.0  2.0  3.5;
                        0.0  2.5  4.0;
                        0.0  3.0  4.5]
    @test output[3] == [0.0  0.0  0.0;
                        0.0  0.0  2.0;
                        0.0  0.0  2.25]
end

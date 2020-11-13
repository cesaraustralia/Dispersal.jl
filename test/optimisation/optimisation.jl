using DynamicGrids, Dispersal, Test, Colors, Flatten, LossFunctions, Optim
using Dispersal: stepfromframe, SingleCoreReplicates, ThreadedReplicates
using DynamicGrids: grid2image

@testset "step from frame" begin
    @test stepfromframe(1, 1, 1) == 1
    @test stepfromframe(1, 5, 1) == 5
    @test stepfromframe(2, 1, 1) == 1
    @test stepfromframe(2, 2, 1) == 1
    @test stepfromframe(2, 3, 1) == 2
    @test stepfromframe(2, 4, 1) == 2
    @test stepfromframe(2, 5, 1) == 3
    @test stepfromframe(3, 4, 1) == 2
    @test stepfromframe(3, 5, 1) == 2
end


@testset "SimpleObjective" begin

    init = [1.0 0.0 1.0
            0.0 0.0 0.0
            1.0 1.0 0.0]

    target = [8.0 0.0 8.0
              0.0 0.0 0.0
              8.0 8.0 0.0]

    # Test the last simulation frame matches the target when the correct rate is used
    tspan = 1:4
    ngroups = 1
    groupsize = 1
    rate = log(2.0)
    ruleset = Ruleset(ExponentialGrowth(intrinsicrate=rate))
    output = ArrayOutput(init; tspan=tspan)
    sim!(output, ruleset)

    @test output[end] == target

    # Test optim finds the correct intrinsic growth rate
    output = ArrayOutput(init; tspan=tspan)
    objective = SimpleObjective(target)
    loss = LogitDistLoss()

    p = Parametriser(ruleset, output, objective, identity, loss, ngroups, groupsize);
    @test DynamicGrids.ruleset(p::Parametriser) === ruleset
    @test Dispersal.output(p::Parametriser) === output
    @test Dispersal.transform(p::Parametriser) == p.transform
    @test Dispersal.objective(p::Parametriser) === objective
    @test Dispersal.loss(p::Parametriser) === loss
    @test Dispersal.ngroups(p::Parametriser) == 1
    @test Dispersal.groupsize(p::Parametriser) == 1
    @test Dispersal.threading(p::Parametriser) == SingleCoreReplicates()
    @test DynamicGrids.tspan(p::Parametriser) == 1:4
    @test Dispersal.data(p::Parametriser) == p.data
    @test Dispersal.targetbuffer(p::Parametriser) == target
    @test Dispersal.targetbuffer(p::Parametriser) !== target
    @test Dispersal.predictionbuffer(p::Parametriser) == target
    @test Dispersal.predictionbuffer(p::Parametriser) !== target

    res = Optim.optimize(p, [log(1.8)], [log(2.2)], [log(1.85)], SAMIN(), Optim.Options(iterations=1000))
    @test res.minimizer[1] ≈ rate atol=4

end

@testset "RegionObjective" begin
    init =  [1.0  1.0  0.5
             0.0  0.0  0.0
             0.0  0.0  0.0]

    regionlookup = [1  2  3
                    1  2  3
                    2  2  0]

    occurance = Bool[0  1  0  # 1
                     1  0  1  # 2
                     1  1  0] # 3

    ngroups = 1
    groupsize = 1
    framesperstep = 2
    steps = 3
    tspan = 1:framesperstep * steps
    detectionthreshold = 0.1
    transform = x -> 2x - 1

    objective = RegionObjective(detectionthreshold, regionlookup, occurance, framesperstep, 1)
    output = RegionOutput(init; tspan=tspan, objective=objective)
    ruleset = Ruleset(ExponentialGrowth(intrinsicrate = log(2.0)))
    loss = ZeroOneLoss()

    @testset "SingleCoreReplicates" begin
        threading = SingleCoreReplicates()
        parametriser = Parametriser(ruleset, output, objective, transform, loss, ngroups, groupsize, threading)
        parametriser(flatten(ruleset))
        res = Optim.optimize(parametriser, [log(1.8)], [log(2.2)], [log(1.85)], SAMIN(), Optim.Options(iterations=1000))
    end

    @testset "ThreadedReplicates" begin
        threading = ThreadedReplicates()
        parametriser = Parametriser(ruleset, output, objective, transform, loss, ngroups, groupsize, threading)
        parametriser(flatten(ruleset))
        res = Optim.optimize(parametriser, [log(1.8)], [log(2.2)], [log(1.85)], SAMIN(), Optim.Options(iterations=1000))
    end

    # TODO test the result with a real dispersal function and predict parameters

end

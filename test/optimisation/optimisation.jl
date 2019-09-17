using DynamicGrids, Dispersal, Test, Colors, Flatten, LossFunctions, Optim
using Dispersal: stepfromframe
using DynamicGrids: frametoimage

@testset "step from frame" begin
    # stepfromframe(framesperstep, frame)
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


# @testset "SimpleObjective" begin

init =  [1.0 0.0 1.0
         0.0 0.0 0.0
         1.0 1.0 0.0]

target = [8.0 0.0 8.0
          0.0 0.0 0.0
          8.0 8.0 0.0]

# Test the last simulation frame matches the target when the correct rate is used
tstop = 4
ngroups = 1
groupsize = 1
rate = log(2.0)
ruleset = Ruleset(ExactExponentialGrowth(intrinsicrate = rate); init=init)
output = ArrayOutput(init, tstop)
sim!(output, ruleset; tstop=tstop)

@test output[end] == target

# Test optim finds the correct intrinsic growth rate
output = ArrayOutput(init, tstop)
objective = SimpleObjective(target)
loss = LogitDistLoss()
parametriser = Parametriser(ruleset, output, objective, identity, loss, ngroups, groupsize, tstop)
res = Optim.optimize(parametriser, [log(1.8)], [log(2.2)], [log(1.85)], SAMIN(), Optim.Options(iterations=1000))
@test res.minimizer[1] ≈ rate atol=4

# @testset "RegionObjective" begin
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
tstop = framesperstep * steps
detectionthreshold = 0.1
transform = x -> 2x - 1


# TODO do this with a real dispersal function and predict parameters
objective = RegionObjective(detectionthreshold, regionlookup, occurance, framesperstep, 1)
output = RegionOutput(init, tstop, objective)
ruleset = Ruleset(ExactExponentialGrowth(intrinsicrate = log(2.0)); init=init)
loss = ZeroOneLoss()
parametriser = Parametriser(ruleset, output, objective, transform, loss, ngroups, groupsize, tstop)

parametriser(flatten(ruleset))

using CellularAutomataBase, Dispersal, Test, Colors, Flatten, LossFunctions
using Dispersal: stepfromframe 
using CellularAutomataBase: frametoimage

@testset "step from frame" begin
    # stepfromframe(framesperstep, frame)
    @test stepfromframe(1, 1) == 1
    @test stepfromframe(1, 5) == 5
    @test stepfromframe(2, 1) == 1
    @test stepfromframe(2, 2) == 1
    @test stepfromframe(2, 3) == 2
    @test stepfromframe(2, 4) == 2
    @test stepfromframe(2, 5) == 3
    @test stepfromframe(3, 4) == 2
    @test stepfromframe(3, 5) == 2
end

@testset "image processor colors regions by fit" begin

    init =  [1.0  1.0  0.5 
             0.0  0.0  0.0 
             0.0  0.0  0.0]

    regionlookup = [1  2  3 
                    1  2  3 
                    2  2  0]

    occurance = Bool[0  1  0  # 1
                     1  0  1  # 2
                     1  1  0] # 3

    image = [RGB24(1.0,0.0,0.0)  RGB24(1.0,1.0,1.0)  RGB24(0.502,0.502,0.502) 
             RGB24(1.0,0.0,1.0)  RGB24(1.0,1.0,0.0)  RGB24(1.0,1.0,0.0) 
             RGB24(1.0,1.0,0.0)  RGB24(1.0,1.0,0.0)  RGB24(0.0,1.0,1.0)]

    framesperstep = 2
    truepostivecolor =   (1.0, 1.0, 1.0)
    falsepositivecolor = (1.0, 0.0, 0.0)
    truenegativecolor =  (1.0, 0.0, 1.0)
    falsenegativecolor = (1.0, 1.0, 0.0)
    maskcolor =          (0.0, 1.0, 1.0)
    objective = RegionObjective(0.0, regionlookup, occurance, framesperstep)

    processor = ColorRegionFit(objective, truepostivecolor, falsepositivecolor,
                               truenegativecolor, falsenegativecolor, maskcolor)

    output = ArrayOutput(init, 1)

    @test image == frametoimage(processor, output, init, 1) 

end

struct TestObjective <: AbstractObjective end

Dispersal.outputtoprediction(obj::TestObjective, output) = 
    output.frames[end]

Dispersal.target(obj::TestObjective) = [1 2 3 
                                        4 5 6
                                        7 8 9]


@testset "Parametriser" begin

    init =  [1.0  1.0  0.5 
             0.0  0.0  0.0 
             0.0  0.0  0.0]

    regionlookup = [1  2  3 
                    1  2  3 
                    2  2  0]

    occurance = Bool[0  1  0  # 1
                     1  0  1  # 2
                     1  1  0] # 3

    nreplicates = 2
    steps = 3
    rule = Ruleset(ExactExponentialGrowth(intrinsicrate = log(2.0)); init=init)
    objective = TestObjective()
    loss = ZeroOneLoss()

    p = Parametriser(rule, objective, loss, nreplicates, 10)
    p(flatten(rule))

    # TODO: actually test the output somehow.
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

    framesperstep = 2
    nreplicates = 2
    detectionthreshold = 0.1
    steps = 3
    rule = Ruleset(ExactExponentialGrowth(intrinsicrate = log(2.0)); init=init)
    objective = RegionObjective(detectionthreshold, regionlookup, occurance, framesperstep)
    loss = ZeroOneLoss()
    p = Parametriser(rule, objective, loss, nreplicates, 6)

    p(flatten(rule))
end

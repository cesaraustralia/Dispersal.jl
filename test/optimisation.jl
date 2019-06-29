using CellularAutomataBase, Dispersal, Test, Colors, Flatten, LossFunctions
using Dispersal: stepfromframe 
using CellularAutomataBase: process_frame

@testset "step from frame" begin
    # stepfromframe(frames_per_step, frame)
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

    region_lookup = [1  2  3 
                     1  2  3 
                     2  2  0]

    occurance = Bool[0  1  0  # 1
                     1  0  1  # 2
                     1  1  0] # 3

    image1 = [RGB24(1.0,0.0,0.0)  RGB24(1.0,1.0,1.0)  RGB24(0.502,0.502,0.502) 
              RGB24(1.0,0.0,1.0)  RGB24(1.0,1.0,0.0)  RGB24(1.0,1.0,0.0) 
              RGB24(1.0,1.0,0.0)  RGB24(1.0,1.0,0.0)  RGB24(0.0,1.0,1.0)]

    image2 = [RGB24(1.0,1.0,1.0)  RGB24(1.0,0.0,0.0)  RGB24(0.502,0.502,0.502) 
              RGB24(1.0,1.0,0.0)  RGB24(1.0,0.0,1.0)  RGB24(1.0,1.0,0.0) 
              RGB24(1.0,0.0,1.0)  RGB24(1.0,0.0,1.0)  RGB24(0.0,1.0,1.0)]

    image3 = [RGB24(1.0,0.0,0.0)  RGB24(1.0,1.0,1.0)  RGB24(0.502,0.0,0.0) 
              RGB24(1.0,0.0,1.0)  RGB24(1.0,1.0,0.0)  RGB24(1.0,0.0,1.0) 
              RGB24(1.0,1.0,0.0)  RGB24(1.0,1.0,0.0)  RGB24(0.0,1.0,1.0)]

    frames_per_step = 2
    truepostivecolor =   (1.0, 1.0, 1.0)
    falsepositivecolor = (1.0, 0.0, 0.0)
    truenegativecolor =  (1.0, 0.0, 1.0)
    falsenegativecolor = (1.0, 1.0, 0.0)
    maskcolor =          (0.0, 1.0, 1.0)

    processor = ColorRegionFit(RegionObjective(0.0, occurance, region_lookup, frames_per_step),
                               truepostivecolor, falsepositivecolor,
                               truenegativecolor, falsenegativecolor, maskcolor)
    # TODO how to test this without GTK?
    # output = GtkOutput(init, processor=processor; fps=1, store=true)

    # @test image1 == process_frame(output, init, 1) 
    # @test image1 == process_frame(output, init, 2) 
    # @test image2 == process_frame(output, init, 3) 
    # @test image2 == process_frame(output, init, 4) 
    # @test image3 == process_frame(output, init, 5) 
    # @test image3 == process_frame(output, init, 6) 

end

struct TestObjective <: AbstractObjective end

Dispersal.outputtoprediction(obj::TestObjective, output) = begin
    output.frames[end]
end
Dispersal.targets(obj::TestObjective) = [1 2 3 
                                         4 5 6
                                         7 8 9]


@testset "Parametriser" begin

    init =  [1.0  1.0  0.5 
             0.0  0.0  0.0 
             0.0  0.0  0.0]

    region_lookup = [1  2  3 
                     1  2  3 
                     2  2  0]

    occurance = Bool[0  1  0  # 1
                     1  0  1  # 2
                     1  1  0] # 3

    num_replicates = 2
    steps = 3
    rule = Ruleset(ExactExponentialGrowth(intrinsicrate = log(2.0)); init=init)
    objective = TestObjective()
    loss = ZeroOneLoss()

    p = Parametriser(rule, objective, loss, num_replicates, 5)
    p(flatten(rule))

    # TODO: actually test the output somehow.
end

@testset "RegionObjective" begin
    init =  [1.0  1.0  0.5 
             0.0  0.0  0.0 
             0.0  0.0  0.0]

    region_lookup = [1  2  3 
                     1  2  3 
                     2  2  0]

    occurance = Bool[0  1  0  # 1
                     1  0  1  # 2
                     1  1  0] # 3

    frames_per_step = 2
    num_replicates = 2
    detection_threshold = 0.1
    steps = 3
    rule = Ruleset(ExactExponentialGrowth(intrinsicrate = log(2.0)); init=init)
    objective = RegionObjective(detection_threshold, region_lookup, occurance, frames_per_step)
    loss = ZeroOneLoss()
    p = Parametriser(rule, objective, loss, num_replicates, 6)

    p(flatten(rule))
end

using Cellular, Dispersal, Test, Colors, Flatten
using Dispersal: step_from_frame 

@testset "step from frame" begin
    # step_from_frame(frames_per_step, t)
    @test step_from_frame(1, 1) == 1
    @test step_from_frame(1, 5) == 5
    @test step_from_frame(2, 1) == 1
    @test step_from_frame(2, 2) == 1
    @test step_from_frame(2, 3) == 2
    @test step_from_frame(2, 4) == 2
    @test step_from_frame(2, 5) == 3
    @test step_from_frame(3, 4) == 2
    @test step_from_frame(3, 5) == 2
end

@testset "SumOutput sums frames" begin
    init =  [1.0  4.0  7.0;
             2.0  5.0  8.0;
             3.0  6.0  9.0]

    test1 = [3.0  12.0  21.0;
             6.0  15.0  24.0;
             9.0  18.0  27.0]

    test2 = [12.0  48.0  84.0;
             24.0  60.0  96.0;
             36.0  72.0 108.0]

    steps = 2
    frames_per_step = 2
    tstop = 4

    model = Models(ExactExponentialGrowth(intrinsicrate = log(2.0)))
    output = SumOutput(init, frames_per_step, steps, Cellular.NullOutput())
    sim!(output, model, init; tstop=tstop)

    @test output[1] == test1
    @test output[2] == test2
end

@testset "image processor colors regions by fit" begin

    init =  [1.0  1.0  0.5;
             0.0  0.0  0.0;
             0.0  0.0  0.0]

    region_lookup = [1  2  3;
                     1  2  3;
                     2  2  0]

    occurance = Bool[0  1  0; # 1
                     1  0  1; # 2
                     1  1  0] # 3

    image1 = [RGB24(1.0,0.0,0.0)  RGB24(1.0,1.0,1.0)  RGB24(0.502,0.502,0.502);
              RGB24(1.0,0.0,1.0)  RGB24(1.0,1.0,0.0)  RGB24(1.0,1.0,0.0);
              RGB24(1.0,1.0,0.0)  RGB24(1.0,1.0,0.0)  RGB24(0.0,1.0,1.0)]

    image2 = [RGB24(1.0,1.0,1.0)  RGB24(1.0,0.0,0.0)  RGB24(0.502,0.502,0.502);
              RGB24(1.0,1.0,0.0)  RGB24(1.0,0.0,1.0)  RGB24(1.0,1.0,0.0);
              RGB24(1.0,0.0,1.0)  RGB24(1.0,0.0,1.0)  RGB24(0.0,1.0,1.0)]

    image3 = [RGB24(1.0,0.0,0.0)  RGB24(1.0,1.0,1.0)  RGB24(0.502,0.0,0.0);
              RGB24(1.0,0.0,1.0)  RGB24(1.0,1.0,0.0)  RGB24(1.0,0.0,1.0);
              RGB24(1.0,1.0,0.0)  RGB24(1.0,1.0,0.0)  RGB24(0.0,1.0,1.0)]

    frames_per_step = 2
    truepostivecolor =   (1.0, 1.0, 1.0)
    falsepositivecolor = (1.0, 0.0, 0.0)
    truenegativecolor =  (1.0, 0.0, 1.0)
    falsenegativecolor = (1.0, 1.0, 0.0)
    maskcolor =          (0.0, 1.0, 1.0)

    processor = ColorRegionFit(frames_per_step, occurance, region_lookup, 
                               truepostivecolor, falsepositivecolor,
                               truenegativecolor, falsenegativecolor, maskcolor)
    output = GtkOutput(init, processor=processor; fps=1, store=true)

    @test image1 == Cellular.process_frame(output, init, 1) 
    @test image1 == Cellular.process_frame(output, init, 2) 
    @test image2 == Cellular.process_frame(output, init, 3) 
    @test image2 == Cellular.process_frame(output, init, 4) 
    @test image3 == Cellular.process_frame(output, init, 5) 
    @test image3 == Cellular.process_frame(output, init, 6) 

end

@testset "RegionParametriser" begin
    init =  [1.0  1.0  0.5;
             0.0  0.0  0.0;
             0.0  0.0  0.0]

    region_lookup = [1  2  3;
                     1  2  3;
                     2  2  0]

    occurance = Bool[0  1  0; # 1
                     1  0  1; # 2
                     1  1  0] # 3

    frames_per_step = 2
    num_replicates = 2
    detection_threshold = 0.1
    steps = 3
    model = Models(ExactExponentialGrowth(intrinsicrate = log(2.0)))
    output = SumOutput(init, frames_per_step, steps, Cellular.NullOutput())

    p = RegionParametriser(output, model, init, occurance, region_lookup, 
                 frames_per_step, num_replicates, detection_threshold)
    p(flatten(model))

    # TODO: actually test the output somehow.
end

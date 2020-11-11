using DynamicGrids, Dispersal, Test, Colors
import DynamicGrids: grid2image

init =  [1.0  1.0  0.5
         0.0  0.0  0.0
         0.0  0.0  0.0]

regionlookup = [1  2  3
                1  2  3
                2  2  0]

occurance = Bool[0  1  0
                 1  0  1
                 1  1  0]

start = 1
detectionthreshold = 0.0
framesperstep = 2
objective = RegionObjective(detectionthreshold, regionlookup, occurance, framesperstep, start)

@testset "region output works" begin
    rule = ExponentialGrowth(intrinsicrate = log(2.0))
    output = RegionOutput(init; objective=objective, tspan=1:5)
    sim!(output, rule)
end

@testset "image processor colors regions by fit" begin
    image = [RGB24(0.2,0.2,0.2)  RGB24(1.0,1.0,1.0)  RGB24(0.5,0.5,0.5)
             RGB24(1.0,0.0,1.0)  RGB24(1.0,1.0,0.0)  RGB24(1.0,1.0,0.0)
             RGB24(1.0,1.0,0.0)  RGB24(1.0,1.0,0.0)  RGB24(0.0,1.0,1.0)]

    truescheme = Greyscale()
    falsescheme = Greyscale(min=nothing, max=0.2)
    truezerocolor =  (1.0, 0.0, 1.0)
    falsezerocolor = (1.0, 1.0, 0.0)
    maskcolor =      (0.0, 1.0, 1.0)

    processor = ColorRegionFit(objective, truescheme, falsescheme,
                               truezerocolor, falsezerocolor, maskcolor)
    output = DynamicGrids.NoDisplayImageOutput(init; tspan=1:5)

    # This output needs an update anyway
    @test image == grid2image(processor, output, Ruleset(), (_default_=init,), 1)
end


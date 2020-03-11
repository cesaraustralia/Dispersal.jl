using Dispersal, Test, DimensionalData
using DimensionalData: DimensionalArray, X, Y, Time

init =  [0.0 0.0 0.0 0.0   0.0 0.0 0.0;
         0.0 0.0 0.0 0.0   0.0 0.0 0.0;
         0.0 0.0 0.0 0.0 100.0 0.0 0.0;
         0.0 0.0 0.0 0.0   0.0 0.0 0.0;
         0.0 0.0 0.0 0.0   0.0 0.0 0.0;
         0.0 0.0 0.0 0.0   0.0 0.0 0.0;
         0.0 0.0 0.0 0.0   0.0 0.0 0.0]

test2 = [0.0  0.0  4.0  4.0  4.0  4.0  4.0;
         0.0  0.0  4.0  4.0  4.0  4.0  4.0;
         0.0  0.0  4.0  4.0  4.0  4.0  4.0;
         0.0  0.0  4.0  4.0  4.0  4.0  4.0;
         0.0  0.0  4.0  4.0  4.0  4.0  4.0;
         0.0  0.0  0.0  0.0  0.0  0.0  0.0;
         0.0  0.0  0.0  0.0  0.0  0.0  0.0]

test3 = [0.48  0.96  1.44  1.92  2.4  1.92  1.44;
         0.64  1.28  1.92  2.56  3.2  2.56  1.92;
         0.8   1.6   2.4   3.2   4.0  3.2   2.4 ;
         0.64  1.28  1.92  2.56  3.2  2.56  1.92;
         0.48  0.96  1.44  1.92  2.4  1.92  1.44;
         0.32  0.64  0.96  1.28  1.6  1.28  0.96;
         0.16  0.32  0.48  0.64  0.8  0.64  0.48]

struct TestFormulation <: KernelFormulation end
(f::TestFormulation)(d) = 1.0

# Dispersal in radius 2 neighborhood
hood = DispersalKernel{2}(; formulation=TestFormulation())

ruleset = Ruleset(InwardsPopulationDispersal(;neighborhood=hood); init=init)
output = ArrayOutput(init, 3)
sim!(output, ruleset; tspan=(1, 3))
@test output[1] == init
@test output[2] == test2
@test output[3] ≈ test3

ruleset = Ruleset(OutwardsPopulationDispersal(;neighborhood=hood); init=init)
output = ArrayOutput(init, 3)
sim!(output, ruleset; tspan=(1, 3))
@test output[1] == init
@test output[2] == test2
@test_broken output[3] ≈ test3


layerdata = cat([1.0  1.0  1.0  1.0  1.0  1.0  1.0
                 1.0  1.0  1.0  1.0  1.0  1.0  1.0
                 1.0  1.0  1.0  1.0  1.0  1.0  1.0
                 1.0  1.0  1.0  1.0  1.0  1.0  1.0
                 1.0  1.0  1.0  1.0  1.0  1.0  1.0
                 1.0  1.0  1.0  1.0  1.0  1.0  1.0
                 1.0  1.0  1.0  1.0  1.0  1.0  1.0],
                [0.0  0.0  0.0  0.0  1.0  1.0  1.0
                 0.0  0.0  0.0  0.0  1.0  1.0  1.0
                 0.0  0.0  0.0  0.0  1.0  1.0  1.0
                 0.0  0.0  0.0  0.0  1.0  1.0  1.0
                 0.0  0.0  0.0  0.0  1.0  1.0  1.0
                 0.0  0.0  0.0  0.0  1.0  1.0  1.0
                 0.0  0.0  0.0  0.0  1.0  1.0  1.0],
                [1.0  1.0  1.0  1.0  0.0  0.0  0.0
                 1.0  1.0  1.0  1.0  0.0  0.0  0.0
                 1.0  1.0  1.0  1.0  0.0  0.0  0.0
                 1.0  1.0  1.0  1.0  0.0  0.0  0.0
                 1.0  1.0  1.0  1.0  0.0  0.0  0.0
                 1.0  1.0  1.0  1.0  0.0  0.0  0.0
                 1.0  1.0  1.0  1.0  0.0  0.0  0.0]; dims=3)

layer = DimensionalArray(layerdata, (X(1:7), Y(1:7), Time(1:3; grid=RegularGrid()))) 

init =  [0.0 0.0 0.0 0.0 0.0 0.0 0.0;
         0.0 0.0 0.0 0.0 0.0 0.0 0.0;
         0.0 0.0 0.0 0.0 0.0 0.0 9.0;
         0.0 0.0 0.0 0.0 0.0 0.0 0.0;
         0.0 0.0 0.0 0.0 0.0 0.0 0.0;
         0.0 9.0 0.0 0.0 0.0 0.0 0.0;
         0.0 0.0 0.0 0.0 0.0 0.0 0.0]

test2 = [0.0  0.0  0.0  0.0  0.0  0.0  0.0;
         0.0  0.0  0.0  0.0  0.0  1.0  1.0;
         0.0  0.0  0.0  0.0  0.0  1.0  1.0;
         0.0  0.0  0.0  0.0  0.0  1.0  1.0;
         0.0  0.0  0.0  0.0  0.0  0.0  0.0;
         0.0  9.0  0.0  0.0  0.0  0.0  0.0;
         0.0  0.0  0.0  0.0  0.0  0.0  0.0]

test3 = [0.0  0.0  0.0  0.0  0.0  0.0  0.0;
         0.0  0.0  0.0  0.0  0.0  1.0  1.0;
         0.0  0.0  0.0  0.0  0.0  1.0  1.0;
         0.0  0.0  0.0  0.0  0.0  1.0  1.0;
         1.0  1.0  1.0  0.0  0.0  0.0  0.0;
         1.0  1.0  1.0  0.0  0.0  0.0  0.0;
         1.0  1.0  1.0  0.0  0.0  0.0  0.0]

threshold = 0.5
hood = DispersalKernel{1}(; formulation=TestFormulation())
switched = SwitchedInwardsPopulationDispersal(; neighborhood=hood, 
                                              threshold=threshold, layer=layer)
ruleset = Ruleset(switched; init=init)

output = ArrayOutput(init, 4)
sim!(output, ruleset; tspan=(1, 4))

@test output[1] == init
@test output[2] == test2
@test output[3] == test3

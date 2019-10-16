using DynamicGrids, Dispersal, Test


init =  [0.0 0.0 0.0 0.0   0.0 0.0 0.0;
         0.0 0.0 0.0 0.0   0.0 0.0 0.0;
         0.0 0.0 0.0 0.0 100.0 0.0 0.0;
         0.0 0.0 0.0 0.0   0.0 0.0 0.0;
         0.0 0.0 0.0 0.0   0.0 0.0 0.0;
         0.0 0.0 0.0 0.0   0.0 0.0 0.0;
         0.0 0.0 0.0 0.0   0.0 0.0 0.0]

test1 = [0.0  0.0  0.0  0.0  0.0  0.0  0.0;
         0.0  0.0  0.0  0.0  0.0  0.0  0.0;
         0.0  0.0  0.0  0.0 100.0 0.0  0.0;
         0.0  0.0  0.0  0.0  0.0  0.0  0.0;
         0.0  0.0  0.0  0.0  0.0  0.0  0.0;
         0.0  0.0  0.0  0.0  0.0  0.0  0.0;
         0.0  0.0  0.0  0.0  0.0  0.0  0.0]

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

struct TestFormulation <: AbstractKernelFormulation end
(f::TestFormulation)(d) = 1.0

# Dispersal in radius 2 neighborhood
hood = DispersalKernel{2}(; formulation=TestFormulation())

rules = Ruleset(InwardsPopulationDispersal(hood); init=init)
output = ArrayOutput(init, 3)
sim!(output, rules; tspan=(1, 3))
output[1]
output[2]
test2
output[3]
test3
@test output[1] == test1
@test output[2] == test2
@test output[3] ≈ test3

rules = Ruleset(OutwardsPopulationDispersal(hood); init=init)
output = ArrayOutput(init, 3)
sim!(output, rules; tspan=(1, 3))
@test output[1] == test1
@test output[2] == test2
@test output[3] ≈ test3

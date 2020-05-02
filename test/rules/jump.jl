using Dispersal, Test, Random

Random.seed!(1234)

init = [0 0 0 0 0;
        0 0 0 0 0;
        0 0 1 0 0;
        0 0 0 0 0;
        0 0 0 0 0]

rules = Ruleset(JumpDispersal(prob_threshold=1.0, spotrange=1); init=init)
output = ArrayOutput(init, 100)
sim!(output, rules; tspan=(1, 100))

output[2] =
    [0 0 0 0 0;
     0 0 0 0 0;
     0 0 1 0 0;
     0 0 0 1 0;
     0 0 0 0 0]

output[3] =
    [0 0 0 0 0;
     0 0 0 0 0;
     0 0 1 0 0;
     0 0 1 1 0;
     0 0 1 0 0]

output[end] =
    [1 1 1 1 1;
     1 1 1 1 1;
     1 1 1 1 1;
     1 1 1 1 1;
     1 1 1 1 1]

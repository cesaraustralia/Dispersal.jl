using Dispersal, Test

init = [1.0 4.0 7.0;
        2.0 5.0 8.0;
        3.0 6.0 9.0]

output = ArrayOutput(init; tspan=1:2)
rule = AlleeExtinction(minfounders = 8.0)
sim!(output, rule)

@test output[1] == [1.0 4.0 7.0;
                    2.0 5.0 8.0;
                    3.0 6.0 9.0]
@test output[2] == [0.0 0.0 0.0;
                    0.0 0.0 8.0;
                    0.0 0.0 9.0]

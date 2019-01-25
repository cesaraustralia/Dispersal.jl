using Cellular, Dispersal, Test

init = [1.0 4.0 7.0;
        2.0 5.0 8.0;
        3.0 6.0 9.0]

mask = [0.0 1.0 0.0;
        0.0 1.0 1.0;
        1.0 1.0 0.0]

output = ArrayOutput(init, 2)

model = Models(Mask(layers=mask))

sim!(output, model, init; tstop=2)

@test output[1] == [1.0 4.0 7.0;
                    2.0 5.0 8.0;
                    3.0 6.0 9.0]

@test output[2] == [0.0 4.0 0.0;
                    0.0 5.0 8.0;
                    3.0 6.0 0.0]

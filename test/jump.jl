using Cellular, Dispersal, Test

init = [0 0 0 0 0;
        0 0 0 0 0;
        0 0 1 0 0;
        0 0 0 0 0;
        0 0 0 0 0]

model = Models(JumpDispersal(prob_threshold=0.0, spotrange=3))
output = ArrayOutput(init, 20)
sim!(output, model, init; tstop=20)

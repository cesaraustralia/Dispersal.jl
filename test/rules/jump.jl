using DynamicGrids, Dispersal, Test

init = [0 0 0 0 0;
        0 0 0 0 0;
        0 0 1 0 0;
        0 0 0 0 0;
        0 0 0 0 0]

rules = Ruleset(JumpDispersal(prob_threshold=0.0, spotrange=3); init=init)
output = ArrayOutput(init, 20)
sim!(output, rules; tspan=(1, 20))

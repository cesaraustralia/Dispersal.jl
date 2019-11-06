using DynamicGrids, Dispersal, Test

sites = BitArray([1 0 0 0 0;
                  0 0 0 0 0;
                  0 0 0 1 0;
                  0 0 0 0 0;
                  0 0 0 0 0])

quarantine = BitArray([0 0 0 0 0;
                       0 0 0 0 0;
                       0 0 0 0 0;
                       0 0 0 0 0;
                       0 0 0 0 0])

pop = [0.0 1.0 0.0 1.0 1.0;
       0.3 1.0 1.0 0.3 1.0;
       1.0 1.0 1.0 4.0 0.1;
       0.0 1.0 0.0 1.5 2.3;
       1.0 1.0 1.0 2.2 3.0]

init = (pop=pop, quarantine=quarantine)

struct AddOne <: CellRule end
DynamicGrids.applyrule(rule::AddOne, data, state, index) = state + one(state)


detection = Detection(population=:pop, quarantine=:quarantine,
                      sites=sites,
                      detection=Dispersal.ThresholdDetection(7),
                      neighborhood=RadialNeighborhood{1}())


ruleset = MultiRuleset(rulesets=(pop=Ruleset(AddOne()), quarantine=Ruleset()),
                    interactions=(detection,),
                    init=init)

output = ArrayOutput(init, 10)
# using DynamicGridsGtk
# processor=DynamicGrids.ThreeColor(colors=(DynamicGrids.Green(), DynamicGrids.Red()))
# output = GtkOutput(init; processor=processor, minval=(0, 0), maxval=(10, 1), store=true)
# output.running = false
sim!(output, ruleset; tspan=(1, 10), fps=3);

@test output[1][:quarantine] == output[3][:quarantine] == [0 0 0 0 0;
                                                           0 0 0 0 0;
                                                           0 0 0 0 0;
                                                           0 0 0 0 0;
                                                           0 0 0 0 0]
# First quarantine
@test output[4][:quarantine] == output[7][:quarantine] == [0 0 0 0 0;
                                                           0 0 1 1 1;
                                                           0 0 1 1 1;
                                                           0 0 1 1 1;
                                                           0 0 0 0 0]
# Second quarantine
@test output[8][:quarantine] == output[10][:quarantine] == [1 1 0 0 0;
                                                            1 1 1 1 1;
                                                            0 0 1 1 1;
                                                            0 0 1 1 1;
                                                            0 0 0 0 0]

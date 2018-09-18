using Revise,
      Dispersal,
      Cellular,
      CuArrays, 
      CUDAnative,
      Test

using Dispersal: hudgins_precalc

@testset "Hudgins" begin

    global init =  Float32[0.0 0.0 1.0 0.0 0.0;
                           0.0 0.0 1.0 0.0 0.0;
                           0.0 0.0 1.0 0.0 0.0;
                           0.0 0.0 0.0 0.0 0.0;
                           0.0 0.0 0.0 0.0 0.0]
    global suit =  Float32[1.0 0.0 1.0 1.0 0.0;
                           0.0 0.0 1.0 1.0 1.0;
                           1.0 1.0 1.0 1.0 0.0;
                           1.0 1.0 0.0 1.0 1.0;
                           1.0 0.0 1.0 1.0 1.0]
    global human = Float32[1.0 1.0 0.0 0.0 0.0;
                           1.0 0.0 0.0 0.0 0.0;
                           0.0 0.0 1.0 0.0 0.0;
                           0.0 0.0 0.0 0.0 1.0;
                           0.0 0.0 0.0 0.0 0.0]

    human = CuArray(human)
    suit = CuArray(suit)
    suit_layer = SuitabilityLayer(suit);
    human_layer = HumanLayer(human);

    precalc = hudgins_precalc(init, suit_layer, human_layer)

    model = Models(HudginsDispersal())
    output = ArrayOutput(init)
    layers = (human_layer, suit_layer)
    sim!(output, model, init, layers, precalc; time=30);
    @test length(output) == 30
end

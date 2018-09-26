using Revise,
      Dispersal,
      Cellular,
      Test
using Dispersal: precalc_human_dispersal, CellMagnitude, CellInterval, build_cell_pop_index

setup(x) = x


@testset "Human dispersal" begin

    global human = setup([0.0 1.0 2.0 3.0 4.0;
                          0.1 1.1 2.1 3.1 4.1;
                          0.2 1.2 2.2 3.2 4.2;
                          0.3 1.3 2.3 3.3 4.3;
                          0.4 1.4 2.4 3.4 4.4])

    global cellsize = 1
    global take = 200
    # human = (rand(400, 400) .* 3) .^ 6

    @time precalc, props = precalc_human_dispersal(human, exponential, (0.05), cellsize, take)
    precalc[100,10]

    Profile.clear()
    @profile 
    using Profile
    using ProfileView
    ProfileView.view()
    precalc[3,5]

    # global layers = SuitabilityLayer(suit)
    # srand(1234)
    # global model = Models(HumanDispersal(prob_threshold=0.5))
    # global output = ArrayOutput(init)
    # sim!(output, model, init, layers; time=3)

    a = [1 2; 3 4]
    b = [1,2,3,4]
    b .= vec(a)
    broadcast!(b -> b, a, b)

end


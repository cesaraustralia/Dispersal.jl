using Dispersal: precalc_human_dispersal, CellMagnitude, CellInterval, build_cell_pop_index, populate!
using Random

@testset "Human dispersal" begin
    Random.seed!(1234)
    global human = setup(ones(5, 5))
    global init = zero(human)
    init[3, 3] = 1.0
    global cellsize = 1
    global take = 9
    global human_exponent = 1
    global dist_exponent = 1
    global precalc, prop = precalc_human_dispersal(human, cellsize, take,
                                                   human_exponent, dist_exponent)
    @test length(precalc) == 25
    @test length(precalc[1, 1]) == take

    a = zeros(5, 5)
    b = zeros(5, 5)

    populate!(a, precalc[1, 1])
    populate!(b, precalc[5, 5])
    @test reverse(a, dims=1) == reverse(b, dims=2)

    a = zeros(Int, 5, 5)
    b = zeros(Int, 5, 5)
    c = zeros(Int, 5, 5)
    d = zeros(Int, 5, 5)
    populate!(a, precalc[1, 1])
    populate!(b, precalc[5, 5])
    populate!(c, precalc[5, 3])
    populate!(d, precalc[3, 3])

    @test a == [1 1 1 0 0;
                1 1 1 0 0;
                1 1 1 0 0;
                0 0 0 0 0;
                0 0 0 0 0]

    @test b == [0 0 0 0 0;
                0 0 0 0 0;
                0 0 1 1 1;
                0 0 1 1 1;
                0 0 1 1 1]

    @test c == [0 0 0 0 0;
                0 0 0 0 0;
                0 0 1 0 0;
                0 1 1 1 0;
                1 1 1 1 1]

    @test d == [0 0 0 0 0;
                0 1 1 1 0;
                0 1 1 1 0;
                0 1 1 1 0;
                0 0 0 0 0]


    global model = Models(HumanDispersal(prob_threshold=0.5, precalc=precalc))
    global output = ArrayOutput(init, 3)
    sim!(output, model, init; tstop=3)
    output[3]
end

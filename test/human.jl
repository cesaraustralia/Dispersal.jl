using Dispersal: precalc_human_dispersal, populate!, MeanDownsampling, upsample_index
using Random

@testset "Human dispersal" begin
    global human_pop = Union{Float32,Missing}[1 for i = 1:5, j = 1:5]
    human_pop[4, 1] = missing
    global init = zeros(Float32, size(human_pop)...)
    init[3, 3] = 10.0
    global cellsize = 1
    global scale = 1
    global aggregator = MeanDownsampling()
    global human_exponent = 1
    global dist_exponent = 1
    global shortlist_len = 9
    global precalc, prop = precalc_human_dispersal(human_pop, cellsize, scale, aggregator,
                        human_exponent, dist_exponent, shortlist_len)
    precalc
    @test length(precalc) == 25
    @test length(precalc[1, 1]) == shortlist_len

    # Propagates missing for missing cells in human pop matrix
    @test ismissing(precalc[4, 1])
    @test ismissing(prop[4, 1])

    a = zeros(5, 5)
    b = zeros(5, 5)

    populate!(a, precalc[1, 1], scale)
    populate!(b, precalc[5, 5], scale)
    @test reverse(a, dims=1) == reverse(b, dims=2)


    a = zeros(Int, 5, 5)
    b = zeros(Int, 5, 5)
    c = zeros(Int, 5, 5)
    d = zeros(Int, 5, 5)
    e = zeros(Int, 5, 5)
    populate!(a, precalc[1, 1], scale)
    populate!(b, precalc[5, 5], scale)
    populate!(c, precalc[5, 3], scale)
    populate!(d, precalc[3, 3], scale)
    populate!(e, precalc[4, 1], scale)

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

    @test e == [0 0 0 0 0;
                0 0 0 0 0;
                0 0 0 0 0;
                0 0 0 0 0;
                0 0 0 0 0]

    global par_a = 1
    global model = Models(HumanDispersal(human_pop, par_a=par_a, cellsize=cellsize, scale=scale, 
                                         aggregator=aggregator, human_exponent=human_exponent, 
                                         dist_exponent=dist_exponent, shortlist_len=shortlist_len
                                        ))
    global output = ArrayOutput(init, 3)

    for i = 1:100
        sim!(output, model, init; tstop=3)
        model.models[1].dispersal_probs[3,3] # probability of a long distance migrant at centre
        single = populate(precalc[3, 3], size(init), scale)
        @test sum(init) == sum(output[3])
        @test sum(init) == sum(output[2])
    end
    # TODO rewrite this it can break with arbitrary code changes 
    @test_broken output[3] ==  [0.0  1.0  0.0  0.0  0.0;
                                0.0  1.0  3.0  0.0  0.0;
                                0.0  0.0  0.0  0.0  0.0;
                                0.0  1.0  0.0  0.0  0.0]
end

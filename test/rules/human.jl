using DynamicGrids, Dispersal, Test, Statistics
using Dispersal: precalc_human_dispersal!, upsample_index, initdownsample,
                 init_dest_shortlist, CellGravity, CellInterval, downsample!

@testset "CellGravity sums and compares gravity" begin
    g = CellGravity(5.0, (1, 1))
    @test g + 3.7 == 8.7
    @test g > 4.0
    @test g < 6.0
    @test g < CellGravity(6.0, (1, 1))
end

@testset "CellGravity sums and compares gravity" begin
    i = CellInterval(5.0, (1, 1))
    @test i > 4.0
    @test i < 6.0
    @test i < CellInterval(6.0, (1, 1))
end


# Setup precalc
human_pop = [1.0f0   2.0f0   0.1f0 0.0f0 1.0f0
             1.0f0   0.0f0   0.3f0 0.4f0 1.0f0
             missing missing 1.0f0 1.0f0 1.0f0
             missing missing 1.0f0 1.0f0 1.0f0
             0.3f0   0.3f0   1.0f0 1.0f0 missing]

test_downsample = [1.0f0   0.2f0 1.0f0
                   missing 1.0f0 1.0f0
                   0.3f0   1.0f0 missing]

human_pop[4, 1] = missing
init = zeros(Float32, size(human_pop)...)
init[3, 3] = 10.0
cellsize = 1
scale = 2
aggregator = mean
human_exponent = 1
dist_exponent = 1
shortlist_len = 4
dispersalperactivity = 0.01
max_dispersers = 4


# Allocate memory
proportion_covered = nothing
human_buffer = initdownsample(human_pop, scale)
dist_buffer = initdownsample(human_pop, scale)
dest_shortlists = init_dest_shortlist(shortlist_len, size(human_buffer))
dist_buffer

precalc_human_dispersal!(dest_shortlists, human_pop, cellsize, scale, aggregator,
                         human_exponent, dist_exponent, shortlist_len, human_buffer, dist_buffer)

@testset "Downsample works" begin
    for i in 1:length(human_buffer)
        d = test_downsample[i]
        h = human_buffer[i]
        ismissing(h) && continue
        @test d == h
    end
end


sl = dest_shortlists[1, 1]
@test length(dest_shortlists) == 9
@test length(dest_shortlists[1, 1]) == 4
@test searchsortedfirst(sl, 1.0) == 4
# With flat human pop a cell mostly disperses to itself
@test sl[4].index == (1, 1)
@test sl[3].index == (2, 2)

# skipmissing(prop) .> 0.0


@testset "Distances are correct" begin
    @test dist_buffer[3, 1] == sqrt(4.0f0^2 + 0)
    @test dist_buffer[3, 2] == sqrt(4.0f0^2 + 2.0f0^2)
end


@testset "precalc is simmetrical with simmetrical inputs" begin
    a = zeros(5, 5)
    b = zeros(5, 5)

    # populate!(a, dest_shortlists[1, 1], scale)
    # populate!(b, dest_shortlists[3, 3], scale)
    # @test reverse(a, dims=1) == reverse(b, dims=2)
end

@testset "human dispersal simulaition maintains total population" begin
    # a = zeros(Int, 5, 5)
    # b = zeros(Int, 5, 5)
    # c = zeros(Int, 5, 5)
    # d = zeros(Int, 5, 5)
    # e = zeros(Int, 5, 5)
    # # populate!(a, precalc[1, 1], scale)
    # # populate!(b, precalc[5, 5], scale)
    # # populate!(c, precalc[5, 3], scale)
    # # populate!(d, precalc[3, 3], scale)
    # # populate!(e, precalc[4, 1], scale)

    # @test a == [1 1 1 0 0;
    #             1 1 1 0 0;
    #             1 1 1 0 0;
    #             0 0 0 0 0;
    #             0 0 0 0 0]

    # @test b == [0 0 0 0 0;
    #             0 0 0 0 0;
    #             0 0 1 1 1;
    #             0 0 1 1 1;
    #             0 0 1 1 1]

    # @test c == [0 0 0 0 0;
    #             0 0 0 0 0;
    #             0 0 1 0 0;
    #             0 1 1 1 0;
    #             1 1 1 1 1]

    # @test d == [0 0 0 0 0;
    #             0 1 1 1 0;
    #             0 1 1 1 0;
    #             0 1 1 1 0;
    #             0 0 0 0 0]

    # @test e == [0 0 0 0 0;
    #             0 0 0 0 0;
    #             0 0 0 0 0;
    #             0 0 0 0 0;
    #             0 0 0 0 0]

    dispersalperpop = 0.1
    rules = Ruleset(HumanDispersal(human_pop, dispersalperpop=dispersalperpop, cellsize=cellsize, scale=scale,
                                  aggregator=aggregator, human_exponent=human_exponent, max_dispersers=max_dispersers,
                                  dist_exponent=dist_exponent, shortlist_len=shortlist_len); init=init)
    output = ArrayOutput(init, 3)

    for i = 1:100
        sim!(output, rules; init=init, tstop=3)
        rules.rules[1].dispersal_probs[3,3] # probability of a long distance migrant at centre
        # Population is allways maintained
        @test sum(init) == sum(output[2])
        @test sum(init) == sum(output[3])
    end

end

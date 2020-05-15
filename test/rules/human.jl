using Dispersal, Test, Statistics
using Dispersal: precalc_human_dispersal!, upsample_index, initdownsample, downsample,
                 CellGravity, CellInterval, downsample!, populate, populate!

@testset "CellGravity sums and compares gravity" begin
    g = CellGravity(5.0, (1, 1))
    @test g + 3.7 == 8.7
    @test g > 4.0
    @test g < 6.0
    @test g < CellGravity(6.0, (1, 1))
end

@testset "CellInterval sums and compares interval" begin
    i = CellInterval(0.5, (1, 1))
    @test i > 0.4
    @test i < 0.6
    @test i < CellInterval(0.6, (1, 1))
end

# Setup precalc
init =      [0.0f0 0.0f0  0.0f0 0.0f0 0.0f0
             0.0f0 0.0f0  0.0f0 0.0f0 0.0f0
             0.0f0 0.0f0 10.0f0 0.0f0 0.0f0
             0.0f0 0.0f0  0.0f0 0.0f0 0.0f0
             0.0f0 0.0f0  0.0f0 0.0f0 0.0f0]
human_pop = [1.0f0   2.0f0   0.1f0 0.0f0 1.0f0
             1.0f0   0.0f0   0.3f0 0.4f0 1.0f0
             missing missing 1.0f0 1.0f0 1.0f0
             missing missing 1.0f0 1.0f0 1.0f0
             0.3f0   0.3f0   1.0f0 1.0f0 missing]

cellsize = 1
scale = 2
aggregator = mean
human_exponent = 2
dist_exponent = 1
nshortlisted = 4
dispersalperactivity = 0.01
max_dispersers = 4

# Allocate memory
human_buffer = downsample(human_pop, aggregator, scale)
distances = initdownsample(human_pop, scale)
gravities, gravity_vector = Dispersal.alloc_gravities(human_buffer, nshortlisted)
dest_shortlists = Dispersal.alloc_dest_shortlist(nshortlisted, size(human_buffer))
gravity_vector

@testset "downsample works" begin
    test_downsample = [1.0f0   0.2f0 1.0f0
                       missing 1.0f0 1.0f0
                       0.3f0   1.0f0 missing]
    for i in 1:length(human_buffer)
        @test human_buffer[i] === test_downsample[i]
    end
    Dispersal.precalc_human_exponent!(human_buffer, human_exponent)
    test_exponent =  [1.0f0    0.2f0^2 1.0f0
                      missing  1.0f0   1.0f0
                      0.09f0   1.0f0   missing]
    for i in 1:length(human_buffer)
        @test human_buffer[i] === test_exponent[i]
    end
end

@testset "distances are correct" begin
    Dispersal.precalc_distances!(distances, dist_exponent, cellsize, scale)
    @test distances[3, 1] == sqrt(2.0f0^2 + 0.0f0^2) * 2
    @test distances[3, 2] == sqrt(2.0f0^2 + 1.0f0^2) * 2
    @test distances[2, 3] == sqrt(1.0f0^2 + 2.0f0^2) * 2
end

@testset "gravities" begin
    Dispersal.precalc_distances!(distances, dist_exponent, cellsize, scale)
    Dispersal.assign_gravities!(gravities, human_buffer, distances, 2, 2)
    @test map(g -> g.gravity, gravities) ==
        [0.35355338f0  0.020000001f0 0.35355338f0
         0.0f0         1.3068552f0     0.5f0
         0.031819806f0 0.5f0         0.0f0]
    @test map(g -> g.index, gravities) ==
        [(1, 1) (1, 2) (1, 3)
         (2, 1) (2, 2) (2, 3)
         (3, 1) (3, 2) (3, 3)]
    gravity_vector .= vec(gravities)
    partialsort!(gravity_vector, nshortlisted, rev=true)
    gravity_shortlist = view(gravity_vector, 1:nshortlisted)
    dest = dest_shortlists[2, 2]
    Dispersal.gravity2inverval!(dest, gravity_shortlist)
    @test map(d -> d.index, dest) == [(1, 1), (2, 3), (3, 2), (2, 2)]
    @test map(d -> d.cumprop, dest) == [0.1328944f0, 0.32083547f0, 0.50877655f0, 1.0f0]
end

@testset "cell" begin
    Dispersal.precalc_cell!(dest_shortlists, nshortlisted, human_buffer, distances, (gravities, gravity_vector), 2, 2)
    dest = dest_shortlists[2, 2]
    @test map(d -> d.index, dest) == [(1, 1), (2, 3), (3, 2), (2, 2)]
    @test map(d -> d.cumprop, dest) == [0.1328944f0, 0.32083547f0, 0.50877655f0, 1.0f0]
end

@testset "shortlists are correctly ordered" begin
    precalc_human_dispersal!(dest_shortlists, human_pop, cellsize, scale, aggregator,
                             human_exponent, dist_exponent, nshortlisted, human_buffer,
                             distances)
    @test length(dest_shortlists) == 9
    @test length(dest_shortlists[1, 1]) == 4
    @test searchsortedfirst(dest_shortlists[1, 1], 1.0) == 4
    @test dest_shortlists[1, 1][4].index == (1, 1)
    @test dest_shortlists[1, 1][3].index == (2, 2)
end

@testset "populate and HumanDispersal construction" begin
    # Setup precalc
    human_pop = Union{Float32,Missing}[1 for i = 1:5, j = 1:5]
    init = zeros(Float32, size(human_pop)...)
    init[3, 3] = 10.0
    cellsize = 1
    scale = 1
    aggregator = mean
    human_exponent = 1
    dist_exponent = 1
    nshortlisted = 9
    dispersalperpop = 0.1

    humandisp = HumanDispersal(human_pop=human_pop, dispersalperpop=dispersalperpop, cellsize=cellsize, scale=scale,
       aggregator=aggregator, human_exponent=human_exponent, max_dispersers=max_dispersers,
       dist_exponent=dist_exponent, nshortlisted=nshortlisted);
    dest_shortlists = humandisp.dest_shortlists

    @testset "precalc is simmetrical with simmetrical inputs" begin
        a = zeros(5, 5)
        b = zeros(5, 5)
        @test (x->x.cumprop).(dest_shortlists[5, 5]) == (x->x.cumprop).(dest_shortlists[1, 1])
        @test (x->x.cumprop).(dest_shortlists[2, 2]) == (x->x.cumprop).(dest_shortlists[4, 4])
        @test a == reverse(reverse(b, dims=2), dims=1)
    end

    a = populate(humandisp, 1, 1)
    b = populate(humandisp, 5, 5)
    c = populate(humandisp[5, 3], size(humandisp), scale)
    d = populate(humandisp[3, 3], size(humandisp), scale)

    @test all(sum.((a, b, c, d)) .≈ 1)

    @test (a .> 0) ==
               [1 1 1 0 0
                1 1 1 0 0
                1 1 1 0 0
                0 0 0 0 0
                0 0 0 0 0]


    @test (b .> 0) ==
               [0 0 0 0 0
                0 0 0 0 0
                0 0 1 1 1
                0 0 1 1 1
                0 0 1 1 1]

    @test (c .> 0) ==
               [0 0 0 0 0
                0 0 0 0 0
                0 0 1 0 0
                0 1 1 1 0
                1 1 1 1 1]

    @test (d .> 0) ==
               [0 0 0 0 0
                0 1 1 1 0
                0 1 1 1 0
                0 1 1 1 0
                0 0 0 0 0]

    @test populate(humandisp, 1, 1) ≈
        reverse(populate(humandisp, 1, 5); dims=2) ≈
        reverse(populate(humandisp, 5, 1); dims=1) ≈
        reverse(reverse(populate(humandisp, 5, 5); dims=1); dims=2)
    @test populate(humandisp, 2, 2) ≈
        reverse(populate(humandisp, 2, 4); dims=2) ≈
        reverse(populate(humandisp, 4, 2); dims=1) ≈
        reverse(reverse(populate(humandisp, 4, 4); dims=1); dims=2)
    @test populate(humandisp, 2, 3) ≈
        rotr90(populate(humandisp, 3, 2)) ≈
        rot180(populate(humandisp, 4, 3)) ≈
        rotl90(populate(humandisp, 3, 4))

    # Symmetry is close but not exactly identical due to truncation
    # of the shortlist and row/column order. This causes small differences
    # choosing the last item in the shortlist when the last index +1
    # has the same gravity. It should make little real difference
    # with larger shortlists as intervals at the end of the shortlist
    # should be small.
    @test (populate(humandisp, 1, 2) ≈
        rotl90(populate(humandisp, 2, 5)) ≈
        rot180(populate(humandisp, 5, 4)) ≈
        rotl90(populate(humandisp, 2, 5))) == false

    combined = populate(humandisp)
    sum(combined) ≈ 25
    @test (combined .> 0) ==
               [1 1 1 1 1
                1 1 1 1 1
                1 1 1 1 1
                1 1 1 1 1
                1 1 1 1 1]

    # Again symmetry is not exact, there are small errors from truncation.
    reverse(combined; dims=1) ≈ combined == false
    reverse(combined; dims=2) ≈ combined == false

    @testset "populate handles missing" begin
        @test populate!(init, missing, scale) === missing
    end
end


@testset "human dispersal simulation maintains total population" begin
    dispersalperpop = 0.1

    humandisp = HumanDispersal(human_pop=human_pop, dispersalperpop=dispersalperpop, cellsize=cellsize, scale=scale,
       aggregator=aggregator, human_exponent=human_exponent, max_dispersers=max_dispersers,
       dist_exponent=dist_exponent, nshortlisted=nshortlisted);

    ruleset = Ruleset(humandisp; init=init)
    output = ArrayOutput(init, 10)

    for i = 1:100
        sim!(output, ruleset; init=init, tspan=(1, 10))
        @test init == output[1]
        @test init != output[10]
        # Population is allways maintained
        @test sum(init) == sum(output[2])
        @test sum(init) == sum(output[3])
        @test sum(init) == sum(output[5])
        @test sum(init) == sum(output[10])
    end

end

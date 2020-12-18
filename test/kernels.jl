using Dispersal, Test, DimensionalData, Setfield, OffsetArrays, DynamicGrids, Random
using DynamicGrids: WritableGridData, SimData, Extent, grids, neighborhood
using Dispersal: dispersalprob

struct TestFormulation <: KernelFormulation end
(f::TestFormulation)(d) = 1.0

# Dispersal in radius 2 neighborhood
hood = DispersalKernel{2}(; formulation=TestFormulation())

@testset "ConstructionBase/Flatten/Setfield compat" begin
    hood2 = @set hood.distancemethod = AreaToArea(10.0)
    @test hood2.distancemethod == AreaToArea(10.0)
    @test DynamicGrids.radius(hood2) == 2
    @test Dispersal.formulation(hood2) == TestFormulation()
    @test hood2.cellsize == 1.0
    @test hood2.kernel == hood.kernel
    @test all(hood2.kernel .== 0.04)
end

@testset "DistanceMethod" begin
    cellsize = 1.0
    f = identity
    x, y = 10, 10

    @test dispersalprob(f, CentroidToCentroid(), x, y, cellsize) == 10sqrt(2)

    subsampling = (1.0, 2.0, 10.0, 20.0, 50.0, 70.0, 80.0)

    proba2c = map(subsampling) do ss
        dispersalprob(f, AreaToCentroid(ss), x, y, cellsize)
    end
    @test proba2c[1] == 10sqrt(2)
    @test issorted(proba2c)
    @test proba2c[end] ≈ proba2c[end-1] atol=1e-5

    proba2a = map(subsampling) do ss
        dispersalprob(f, AreaToArea(ss), x, y, cellsize)
    end
    @test proba2a[1] == 10sqrt(2)
    @test issorted(proba2a)
    @test proba2a[end] ≈ proba2a[end-1] atol=1e-5
end

@testset "buildkernel" begin
    #buildkernel(formulation, distancemethod, cellsize, r)
end

init =  [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
         0.0 0.0 0.0 0.0 81. 0.0 0.0 0.0;
         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]

test2 = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
         0.0 0.0 0.0 9.0 9.0 9.0 0.0 0.0;
         0.0 0.0 0.0 9.0 9.0 9.0 0.0 0.0;
         0.0 0.0 0.0 9.0 9.0 9.0 0.0 0.0;
         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]

test3 = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
         0.0 0.0 1.0 2.0 3.0 2.0 1.0 0.0;
         0.0 0.0 2.0 4.0 6.0 4.0 2.0 0.0;
         0.0 0.0 3.0 6.0 9.0 6.0 3.0 0.0;
         0.0 0.0 2.0 4.0 6.0 4.0 2.0 0.0;
         0.0 0.0 1.0 2.0 3.0 2.0 1.0 0.0;
         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]

@testset "inwards" begin
    hood = DispersalKernel{1}(; formulation=TestFormulation())
    ruleset = Ruleset(InwardsDispersal(;neighborhood=hood))
    output = REPLOutput(init; tspan=1:3, store=true)
    sim!(output, ruleset)
    @test output[1] == init
    @test output[2] == test2
    @test output[3] == test3
end

@testset "outwards" begin
    hood = DispersalKernel{1}(; formulation=TestFormulation())
    ruleset = Ruleset(OutwardsDispersal(;neighborhood=hood))
    output = REPLOutput(init; tspan=1:3, store=true)
    sim!(output, ruleset)
    @test output[1] == init
    @test output[2] == test2
    @test output[3] == test3
end

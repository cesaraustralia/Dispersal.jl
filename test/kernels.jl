using Dispersal
using Test, Setfield, DynamicGrids, StaticArrays
using Dispersal: dispersalprob

struct TestFormulation <: KernelFormulation end
(f::TestFormulation)(d) = 1.0

# Dispersal in radius 2 neighborhood
hood = DispersalKernel(; formulation=TestFormulation())
@test hood == DispersalKernel(Window{1}(), SMatrix{3,3}(fill(1/9, 3, 3)), TestFormulation(), 1.0, CentroidToCentroid())

@testset "ConstructionBase/Flatten/Setfield compat" begin
    hood2 = @set hood.distancemethod = AreaToArea(10.0);
    @test hood2.distancemethod == AreaToArea(10.0)
    @test DynamicGrids.radius(hood2) == 1
    @test Dispersal.formulation(hood2) == TestFormulation()
    @test hood2.cellsize == 1.0
    @test hood2.kernel == hood.kernel
    @test all(hood2.kernel .== 1/9)
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

@testset "Window" begin

    init =  [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
             0.0 0.0 0.0 0.0 81. 0.0 0.0 0.0
             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]

    test2 = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
             0.0 0.0 0.0 9.0 9.0 9.0 0.0 0.0
             0.0 0.0 0.0 9.0 9.0 9.0 0.0 0.0
             0.0 0.0 0.0 9.0 9.0 9.0 0.0 0.0
             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]

    test3 = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
             0.0 0.0 1.0 2.0 3.0 2.0 1.0 0.0
             0.0 0.0 2.0 4.0 6.0 4.0 2.0 0.0
             0.0 0.0 3.0 6.0 9.0 6.0 3.0 0.0
             0.0 0.0 2.0 4.0 6.0 4.0 2.0 0.0
             0.0 0.0 1.0 2.0 3.0 2.0 1.0 0.0
             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]

    @testset "inwards" begin
        ruleset = Ruleset(InwardsDispersal(; radius=1, formulation=TestFormulation()))
        output = REPLOutput(init; tspan=1:3, store=true)
        sim!(output, ruleset)
        @test output[1] == init
        @test output[2] == test2
        @test output[3] == test3
        sim!(output, ruleset; opt=SparseOpt())
        @test output[1] == init
        @test output[2] == test2
        @test output[3] == test3
        sim!(output, ruleset; proc=ThreadedCPU(), opt=SparseOpt())
        @test output[1] == init
        @test output[2] == test2
        @test output[3] == test3
    end

    @testset "outwards" begin
        ruleset = Ruleset(OutwardsDispersal(; radius=1, formulation=TestFormulation()))
        output = REPLOutput(init; tspan=1:3, store=true)
        sim!(output, ruleset)
        @test output[1] == init
        @test output[2] == test2
        @test output[3] == test3
        sim!(output, ruleset; opt=SparseOpt())
        @test output[1] == init
        @test output[2] == test2
        @test output[3] == test3
        sim!(output, ruleset; proc=ThreadedCPU())
        @test output[1] == init
        @test output[2] == test2
        @test output[3] == test3
        sim!(output, ruleset; proc=ThreadedCPU(), opt=SparseOpt())
        @test output[1] == init
        @test output[2] == test2
        @test output[3] == test3
    end

end

@testset "Positional" begin

    init =  [0.0 0.0 0.0 0.0  0.0 0.0 0.0 0.0
             0.0 0.0 0.0 0.0  0.0 0.0 0.0 0.0
             0.0 0.0 0.0 0.0  0.0 0.0 0.0 0.0
             0.0 0.0 0.0 0.0 25.0 0.0 0.0 0.0
             0.0 0.0 0.0 0.0  0.0 0.0 0.0 0.0
             0.0 0.0 0.0 0.0  0.0 0.0 0.0 0.0
             0.0 0.0 0.0 0.0  0.0 0.0 0.0 0.0
             0.0 0.0 0.0 0.0  0.0 0.0 0.0 0.0]

    test2 = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
             0.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0
             0.0 0.0 0.0 5.0 5.0 5.0 0.0 0.0
             0.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0
             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]

    test3 = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
             0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0
             0.0 0.0 0.0 2.0 2.0 2.0 0.0 0.0
             0.0 0.0 1.0 2.0 5.0 2.0 1.0 0.0
             0.0 0.0 0.0 2.0 2.0 2.0 0.0 0.0
             0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0
             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
             0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]

    hood = Positional((-1, 0), (0, -1), (0, 0), (1, 0), (0, 1)) 

    @testset "inwards" begin
        rule = InwardsDispersal(; neighborhood=hood, formulation=TestFormulation())
        output = REPLOutput(init; tspan=1:3, store=true)
        for kw in ((), (opt=SparseOpt(),), (proc=ThreadedCPU(),), (proc=ThreadedCPU(), opt=SparseOpt()))
            sim!(output, rule; kw...)
            @test output[1] == init
            @test output[2] == test2
            @test output[3] == test3
        end
    end

    @testset "outwards" begin
        rule = OutwardsDispersal(; neighborhood=hood, formulation=TestFormulation())
        output = REPLOutput(init; tspan=1:3, store=true)
        for kw in ((), (opt=SparseOpt(),), (proc=ThreadedCPU(),), (proc=ThreadedCPU(), opt=SparseOpt()))
            sim!(output, rule; kw...)
            @test output[1] == init
            @test output[2] == test2
            @test output[3] == test3
        end
    end

end

@testset "kernel formulation" begin
    @testset "ExponentialKernel" begin
        hood = DispersalKernel(; formulation=ExponentialKernel(0.5))
        @test DynamicGrids.radius(hood) == 1
        @test sum(DynamicGrids.kernel(hood)) ≈ 1 atol=10e-10
    end

    @testset "GeometricKernel" begin
        hood = DispersalKernel(; formulation=GeometricKernel(0.5))
        @test DynamicGrids.radius(hood) == 1
        @test sum(DynamicGrids.kernel(hood)) ≈ 1 atol=10e-10
    end

    @testset "GaussianKernel" begin
        hood = DispersalKernel(; formulation=GaussianKernel())
        @test DynamicGrids.radius(hood) == 1
        @test sum(DynamicGrids.kernel(hood)) ≈ 1 atol=10e-10
    end

    @testset "WeibullKernel" begin
        hood = DispersalKernel(; formulation=WeibullKernel())
        @test DynamicGrids.radius(hood) == 1
        @test sum(DynamicGrids.kernel(hood)) ≈ 1 atol=10e-10
    end
end

   

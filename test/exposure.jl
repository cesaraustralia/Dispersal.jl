using Dispersal
using DimensionalData
using Dates
using Test

@testset "Exposure trigger by threshold" begin

    init = (
        exposure = [0.0 0.0 0.0; 
                     0.0 0.0 0.0;
                     0.0 0.0 0.0],
        population = [1000.0 100.0  0.0;
                      0.0    1000.0 0.0;
                      0.0    0.0    10.0],
    )

    active = [0.0 1.0  1.0; 
            0.0 10.0 0.0;
            0.0 0.0  10.0]
    
    exposurerule = ThresholdExposure{Tuple{:exposure,:population},:exposure}(
        active=Aux(:active), pulselevel=20, popthreshold = 100.0, timestep=1
    )
    degradationrule = ExponentialDecrease{:exposure,:exposure}(rate = 0.5, timestep=1)

    @testset "Threshold Exposure" begin        
        # no degradation
        output = ArrayOutput(init; tspan=1:3,  aux=(active=active,))
        sim!(output, exposurerule)
        @test output[1].exposure == [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
        @test output[2].exposure == [0.0 20.0 0.0; 0.0 20.0 0.0; 0.0 0.0 0.0]
        @test output[3].exposure == [0.0 40.0 0.0; 0.0 40.0 0.0; 0.0 0.0 0.0]
        # Degradation
        output = ArrayOutput(init; tspan=1:3,  aux=(active=active,))
        sim!(output, Ruleset(exposurerule, degradationrule))
        @test output[1].exposure == [0.0 0.0     0.0; 0.0 0.0     0.0; 0.0 0.0 0.0]
        @test output[2].exposure ≈  [0.0 12.1306 0.0; 0.0 12.1306 0.0; 0.0 0.0 0.0] atol=0.01
        @test output[3].exposure ≈  [0.0 19.4882 0.0; 0.0 19.4882 0.0; 0.0 0.0 0.0] atol=0.01
        # Degradation
        output = ArrayOutput(init; tspan=1:3,  aux=(active=active,))
        sim!(output, Ruleset(degradationrule, exposurerule))
        @test output[1].exposure == [0.0 0.0     0.0; 0.0 0.0     0.0; 0.0 0.0 0.0]
        @test output[2].exposure == [0.0 20.0    0.0; 0.0 20.0    0.0; 0.0 0.0 0.0]
        @test output[3].exposure ≈  [0.0 32.1306 0.0; 0.0 32.1306 0.0; 0.0 0.0 0.0] atol=0.01
    end

    @testset "variable population" begin
        populationrule = AlleeExtinction{:population,:population}(minfounders = 200.0)
        output = ArrayOutput(init; tspan=1:3,  aux=(active=active,))
        sim!(output, Ruleset(degradationrule, exposurerule , populationrule))
        @test output[1] == (
            exposure = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0],
            population = [1000.0 100.0 0.0; 0.0 1000.0 0.0; 0.0 0.0 10.0])
        @test output[2] == (
            exposure = [0.0 20.0 0.0; 0.0 20.0 0.0; 0.0 0.0 0.0],
            population = [1000.0 0.0 0.0; 0.0 1000.0 0.0; 0.0 0.0 0.0])
        @test output[3] == (
            # only degradation where the population has been removed 
            exposure = [0.0 12.130613194252668 0.0; 0.0 32.13061319425267 0.0; 0.0 0.0 0.0],
            population = [1000.0 0.0 0.0; 0.0 1000.0 0.0; 0.0 0.0 0.0])
    end

    @testset "variable active, constant population" begin
        active = reshape(
            [0.0 1.0 1.0   0.0 0.0 0.0  0.0 1.0 1.0;
            0.0 10.0 0.0  0.0 0.0 0.0  0.0 10.0 0.0;
            0.0 0.0 10.0  0.0 0.0 0.0  0.0 0.0 10.0],
            (3,3,3))
        active_da = DimArray(active, (X, Y, Ti))
        output = ArrayOutput(init; tspan=1:3,  aux=(active=active_da,))
        sim!(output, Ruleset(degradationrule, exposurerule))
        @test output[Ti(1)].exposure == [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
        @test output[Ti(2)].exposure == [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
        @test output[Ti(3)].exposure == [0.0 20.0 0.0; 0.0 20.0 0.0; 0.0 0.0 0.0]
    end
end   

@testset "Exposure following rotation of exposures" begin

    @testset "rotation structure" begin
        rotation(x,y) = Dispersal.Rotation(x,y)

        @test rotation(3,2) + rotation(4,8) == rotation(7,2)
        @test rotation(4,8) + rotation(3,2) == rotation(7,8)

        _updaterotation(x,y) = Dispersal._updaterotation(x, y) 
        @test _updaterotation(rotation(3,2),1) == rotation(4,1)
        @test _updaterotation(rotation(3,2), 4) == rotation(4,4)
        @test _updaterotation(
            rotation(3,DateTime(2020,2,1)), DateTime(2021,1,10)
        ) == rotation(4,DateTime(2021,1,10))

        @test initrotation(1, (3,3)) == [
            rotation(1,1) rotation(1,1) rotation(1,1);
            rotation(1,1) rotation(1,1) rotation(1,1);
            rotation(1,1) rotation(1,1) rotation(1,1)
        ]

        @test initrotation(DateTime(2020,2,1), (2,2)) == [
            rotation(1,DateTime(2020,2,1)) rotation(1,DateTime(2020,2,1)) ;
            rotation(1,DateTime(2020,2,1)) rotation(1,DateTime(2020,2,1))
        ]
    end

    @testset "generate image with rotation structure" begin
        @test DynamicGrids.to_rgb(ObjectScheme(), RS(99.0, 1.0) / 99) == ARGB32(1.0)
        @test DynamicGrids.to_rgb(ObjectScheme(), RS(00.0, 0.0) / 99) == ARGB32(0.0)
        @test DynamicGrids.to_rgb(ObjectScheme(), DynamicGrids.normalise(RS(99.0, 1.0), nothing, 99)) == ARGB32(1.0)
    end
  
    @testset "Rotation rules" begin
        nr, nc = 1, 1
        init = (
            exposure1 = zeros(nr,nc),
            exposure2 = zeros(nr,nc),
            exposure3 = zeros(nr,nc),
            rotation = initrotation(1, (nr,nc)),
            population = ones(nr,nc).*10^4, )

        active = ones(nr,nc)
        
        output = ArrayOutput(init; tspan=1:7,  aux=(active=active,))
        exposurerule1 = RotationExposure{Tuple{:exposure1,:population,:rotation},Tuple{:exposure1,:rotation}}(
            active=Aux(:active), rotationsize=3, rotationindex= 1,
            pulselevel=20.0, popthreshold = 1000.0,  timestep = 1
        )
        exposurerule2 = RotationExposure{Tuple{:exposure2,:population,:rotation},Tuple{:exposure2,:rotation}}(
            active=Aux(:active), rotationsize=3, rotationindex= 2,
            pulselevel=20.0, popthreshold = 1000.0, timestep = 1
        )
        exposurerule3 = RotationExposure{Tuple{:exposure3,:population,:rotation},Tuple{:exposure3,:rotation}}(
            active=Aux(:active), rotationsize=3, rotationindex= 3,
            pulselevel=20.0, popthreshold = 1000.0, timestep = 1
        )
        degradationrule1 = ExponentialDecrease{:exposure1,:exposure1}(rate = 0.5, timestep=1)
        degradationrule2 = ExponentialDecrease{:exposure2,:exposure2}(rate = 0.5, timestep=1)
        degradationrule3 = ExponentialDecrease{:exposure3,:exposure3}(rate = 0.5, timestep=1)

        rule = Ruleset(Chain(
            degradationrule1, exposurerule1,
            degradationrule2, exposurerule2,
            degradationrule3, exposurerule3))
        sim!(output, rule)

        exposure_vals(i) = [output[i].exposure1[1] , output[i].exposure2[1], output[i].exposure3[1]]

        # 1st step: no exposure
        @test exposure_vals(1) == [0.0, 0.0, 0.0]
        # 2nd step: pulse only exposure 1
        @test exposure_vals(2) == [20.0, 0.0, 0.0]
        # 3rd step: degrad exposure 1, pulse exposure 2, 0 exposure 3
        @test exposure_vals(3) ≈ [12.1306, 20.0, 0.0] atol=0.01
        # 4th step:  degrad exposure 1, degrad exposure 2, pulse exposure 3
        @test exposure_vals(4) ≈ [7.3575,  12.1306, 20.0] atol=0.01
        # 5th step:  pulse exposure 1, degrad exposure 2, degrad exposure 3
        @test exposure_vals(5) ≈ [24.4626, 7.3575, 12.1306] atol=0.01
        # 6th step:  degrad exposure 1, pulse exposure 2, degrad exposure 3
        @test exposure_vals(6) ≈ [14.8373, 24.4626, 7.3575] atol=0.01
    end   
end

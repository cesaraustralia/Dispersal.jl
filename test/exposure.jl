using Dispersal
using DynamicGrids, Unitful, Dates
using DimensionalData
using Test

@testset "Test exposure population" begin

    init = (
        pesticide =  [0.0 0.0 0.0; 
                    0.0 0.0 0.0;
                    0.0 0.0 0.0],
        population = [100000.0 100.0 0.0;
                    0.0 1000.0 0.0;
                    0.0 0.0 100000.0], )

    crop = [0.0 1.0 1.0; 
    0.0 10.0 0.0;
    0.0 0.0 10.0]
        
    output = ArrayOutput(init; tspan=1:3,  aux=(crop=crop,))
    ruleTreatment = Threshold_Exposure{Tuple{:pesticide,:population}, :pesticide}(crop=Aux(:crop), pulseLevel=20, popThreshold = 1000)
    rule = Ruleset(ruleTreatment)
    sim!(output, rule)
    @test output[Ti(1)] == (
        pesticide = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0],
        population = [100000.0 100.0 0.0; 0.0 1000.0 0.0; 0.0 0.0 100000.0])
    @test output[Ti(2)] == (
        pesticide = [0.0 0.0 0.0; 0.0 20.0 0.0; 0.0 0.0 20.0],
        population = [100000.0 100.0 0.0; 0.0 1000.0 0.0; 0.0 0.0 100000.0])
    @test output[Ti(3)] == (
        pesticide = [0.0 0.0 0.0; 0.0 32.13061319425267 0.0; 0.0 0.0 32.13061319425267],
        population = [100000.0 100.0 0.0; 0.0 1000.0 0.0; 0.0 0.0 100000.0])

    # If cropping change allong time
    crop = reshape( # TODO: something that change with time
        repeat([0.0 1.0 1.0; 
                0.0 10.0 0.0;
                0.0 0.0 10.0], outer=(1,3) ), (3,3,3))
    cropDD = DimArray(crop, (X, Y, Ti));

    output = ArrayOutput(init; tspan=1:3,  aux=(cropDD=cropDD,))
    
    ruleTreatment = Threshold_Exposure{Tuple{:pesticide,:population}, :pesticide}(crop=Aux(:cropDD), pulseLevel=20, popThreshold = 1000)

    rule = Ruleset(ruleTreatment)
    sim!(output, rule)
    @test output[Ti(1)] == (
        pesticide = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0],
        population = [100000.0 100.0 0.0; 0.0 1000.0 0.0; 0.0 0.0 100000.0])
    @test output[Ti(2)] == (
        pesticide = [0.0 0.0 0.0; 0.0 20.0 0.0; 0.0 0.0 20.0],
        population = [100000.0 100.0 0.0; 0.0 1000.0 0.0; 0.0 0.0 100000.0])
    @test output[Ti(3)] == (
        pesticide = [0.0 0.0 0.0; 0.0 32.13061319425267 0.0; 0.0 0.0 32.13061319425267],
        population = [100000.0 100.0 0.0; 0.0 1000.0 0.0; 0.0 0.0 100000.0])

end   

##################
# GENERIC ROTATION
##################

@testset "rotation structure" begin
    RS(x,y) = Dispersal.RotationStruct(x,y)
   
    @test mod(RS(1,2),2) == mod(1,2)
    @test mod(RS(2,2),2) == mod(2,2)
    @test mod(RS(2,3),2) == mod(2,2)
    @test mod(RS(3,2),2) == mod(3,2)

    @test RS(3,2) + 1 == RS(4,2)
    @test RS(3,2) + 8.3 == RS(11.3,2)

    US(x,y) = Dispersal.updateStep(x, y)
    @test US(RS(3,2), 4) == RS(3,4)
    @test US(RS(3,DateTime(2020,2,1)), DateTime(2021,1,10)) == RS(3,DateTime(2021,1,10))

    @test Rotation_InitGrid(1, (3,3)) == [
        RS(1,1) RS(1,1) RS(1,1);
        RS(1,1) RS(1,1) RS(1,1);
        RS(1,1) RS(1,1) RS(1,1)
    ]

    @test Rotation_InitGrid(DateTime(2020,2,1), (2,2)) == [
        RS(1,DateTime(2020,2,1)) RS(1,DateTime(2020,2,1)) ;
        RS(1,DateTime(2020,2,1)) RS(1,DateTime(2020,2,1))
    ]
end


@testset "Test rotation in exposure" begin
    nr, nc = 1, 1
    init = (
        pesticide1 = zeros(nr,nc),
        pesticide2 = zeros(nr,nc),
        pesticide3 = zeros(nr,nc),
        rotation = Rotation_InitGrid(1, (nr,nc)),
        population = ones(nr,nc).*10^4, )

    crop = ones(nr,nc)
       
    output = ArrayOutput(init; tspan=1:7,  aux=(crop=crop,))
    ruleTreatment1 = Rotation_Exposure{Tuple{:pesticide1,:population,:rotation}, Tuple{:pesticide1,:rotation}}(
        crop=Aux(:crop), rotationSize=3, rotationNumber= 1, pulseLevel=20.0, popThreshold = 1000.0, degradationRate=0.5
    )
    ruleTreatment2 = Rotation_Exposure{Tuple{:pesticide2,:population,:rotation}, Tuple{:pesticide2,:rotation}}(
        crop=Aux(:crop), rotationSize=3, rotationNumber= 2, pulseLevel=20.0, popThreshold = 1000.0,degradationRate=0.5
    )
    ruleTreatment3 = Rotation_Exposure{Tuple{:pesticide3,:population,:rotation}, Tuple{:pesticide3,:rotation}}(
        crop=Aux(:crop), rotationSize=3, rotationNumber= 3, pulseLevel=20.0, popThreshold = 1000.0,degradationRate=0.5
    )
    rule = Ruleset(Chain(ruleTreatment1, ruleTreatment2,ruleTreatment3))
    sim!(output, rule)

    CheckReturn(i) = [output[i].pesticide1[1] , output[i].pesticide2[1], output[i].pesticide3[1]]

    # 1st step: no treatment
    @test CheckReturn(1) == [0.0, 0.0, 0.0]
    # 2nd step: pulse only treatment 1
    @test CheckReturn(2) == [20.0, 0.0, 0.0]
    # 3rd step: degrad treatment 1, pulse treatment 2, 0 treatment 3
    @test CheckReturn(3) ≈ [12.1306, 20.0, 0.0] atol=0.01
    # 4th step:  degrad treatment 1, degrad treatment 2, pulse treatment 3
    @test CheckReturn(4) ≈ [7.3575,  12.1306, 20.0] atol=0.01
    # 5th step:  pulse treatment 1, degrad treatment 2, degrad treatment 3
    @test CheckReturn(5) ≈ [24.4626, 7.3575, 12.1306] atol=0.01
    # 6th step:  degrad treatment 1, pulse treatment 2, degrad treatment 3
    @test CheckReturn(6) ≈ [14.8373, 24.4626, 7.3575] atol=0.01
end   
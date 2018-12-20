@testset "floating inwards dispersal simulation" begin


    global init = setup([  0.0 0.0 0.0 0.0   0.0 0.0 0.0;
                           0.0 0.0 0.0 0.0   0.0 0.0 0.0;
                           0.0 0.0 0.0 0.0 100.0 0.0 0.0;
                           0.0 0.0 0.0 0.0   0.0 0.0 0.0;
                           0.0 0.0 0.0 0.0   0.0 0.0 0.0;
                           0.0 0.0 0.0 0.0   0.0 0.0 0.0;
                           0.0 0.0 0.0 0.0   0.0 0.0 0.0])

    global test1 = setup([ 0.0 0.0 0.0 0.0   0.0 0.0 0.0;
                           0.0 0.0 0.0 0.0   0.0 0.0 0.0;
                           0.0 0.0 0.0 0.0 100.0 0.0 0.0;
                           0.0 0.0 0.0 0.0   0.0 0.0 0.0;
                           0.0 0.0 0.0 0.0   0.0 0.0 0.0;
                           0.0 0.0 0.0 0.0   0.0 0.0 0.0;
                           0.0 0.0 0.0 0.0   0.0 0.0 0.0])

    global test2 = setup([ 0.0  0.0  4.0  4.0  4.0  4.0  4.0;
                           0.0  0.0  4.0  4.0  4.0  4.0  4.0;
                           0.0  0.0  4.0  4.0  4.0  4.0  4.0;
                           0.0  0.0  4.0  4.0  4.0  4.0  4.0;
                           0.0  0.0  4.0  4.0  4.0  4.0  4.0;
                           0.0  0.0  0.0  0.0  0.0  0.0  0.0;
                           0.0  0.0  0.0  0.0  0.0  0.0  0.0])

    global test3 = setup([0.48  0.96  1.44  1.92  2.4  1.92  1.44;
                         0.64  1.28  1.92  2.56  3.2  2.56  1.92;
                         0.8   1.6   2.4   3.2   4.0  3.2   2.4 ;
                         0.64  1.28  1.92  2.56  3.2  2.56  1.92;
                         0.48  0.96  1.44  1.92  2.4  1.92  1.44;
                         0.32  0.64  0.96  1.28  1.6  1.28  0.96;
                         0.16  0.32  0.48  0.64  0.8  0.64  0.48])

    # Dispersal in radius 2 neighborhood

    global hood = DispersalKernel(; f=(d,a)->1.0, radius=2)
    global model = Models(OutwardsPopulationDispersal(neighborhood=hood))
    global output = ArrayOutput(init, 3)
    sim!(output, model, init; tstop=3)
    @test output[1] == test1
    @test output[2] == test2
    @test output[3] ≈ test3

end

using Dispersal, Statistics, Test
using Dispersal: downsample, upsample_index, downsample_index

a = Union{Missing,Float64}[missing 2.0 3.0 4.0;
                           missing 2.0 3.0 4.0;
                           1.0 2.0 3.0 missing;
                           1.0 2.0 3.0 6.0]
scale = 2

@test downsample(a, sum, scale) == [4.0 14.0;
                                    6.0 12.0]

@test downsample(a, mean, scale) == [2.0 3.5;
                                     1.5 4.0]

@test downsample(a, maximum, scale) == [2.0 4.0;
                                        2.0 6.0]

@test downsample(a, minimum, scale) == [2.0 3.0;
                                        1.0 3.0]


@test upsample_index((1, 1), 2) == (1, 1)
@test upsample_index((2, 1), 2) == (3, 1)
@test upsample_index((1, 4), 2) == (1, 7)
@test upsample_index((3, 4), 2) == (5, 7)
@test upsample_index((100, 97), 2) == (199, 193)

@test downsample_index((1, 1), 2) == (1, 1)
@test downsample_index((3, 1), 2) == (2, 1)
@test downsample_index((1, 7), 2) == (1, 4)
@test downsample_index((6, 8), 2) == (3, 4)
@test downsample_index((199, 193), 2) == (100, 97)

@test downsample_index(upsample_index((99,324), 2), 2) == (99,324)



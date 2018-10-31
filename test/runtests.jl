using Revise,
      Cellular,
      Dispersal,
      Test
using Dispersal: suitability, cyclic, sequence_interpolate, neighbors

setup(x) = x

# For manual testing on CUDA
# using CuArrays, CUDAnative
# setup(x) = CuArray(x)

include("growth.jl")
include("dispersal.jl")
include("layers.jl")

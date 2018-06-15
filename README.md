# Dispersal

[![Build Status](https://travis-ci.org/rafaqz/Dispersal.jl.svg?branch=master)](https://travis-ci.org/rafaqz/Dispersal.jl)
[![Coverage Status](https://coveralls.io/repos/rafaqz/Dispersal.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/rafaqz/Dispersal.jl?branch=master)
[![codecov.io](http://codecov.io/github/rafaqz/Dispersal.jl/coverage.svg?branch=master)](http://codecov.io/github/rafaqz/Dispersal.jl?branch=master)


This package builds on [Cellular.jl](https://github.com/rafaqz/Cellular.jl) to
provide a framework for cellular dispersal simulations.

To run a simulation:

```julia
using Cellular
using Dispersal
using Tk

# add a 0.0-1.0 scaled raster array that represents habitat suitability
suit = your_2d_array

# Define a dispersal kernel. This can be any function accepts one argument: 
# distance from the central cell, and returns a value from 0.0 to 1.0.
f = d -> exponential(d, 1)

# Define the neighborhood as our dispersal kernel
hood = DispersalNeighborhood(; f=f, radius=3)

# Define additional raster layers
layers = SuitabilityLayer(suit)

# Define disersal modules
localdisp = LocalDispersal(layers=layers, neighborhood=hood)
jumpdisp = JumpDispersal(layers=layers)

# Build the complete model
model = (localdisp, jumpdisp)

# Create the source array and seed it
source = zeros(Int8, size(suit))
source[24, 354] = 1

# Set the output type
output = TkOutput(source)

# Run the simulation
sim!(source, model, output) 
```

# Dispersal

[![Build Status](https://travis-ci.org/rafaqz/Dispersal.jl.svg?branch=master)](https://travis-ci.org/rafaqz/Dispersal.jl)
[![Coverage Status](https://coveralls.io/repos/rafaqz/Dispersal.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/rafaqz/Dispersal.jl?branch=master)
[![codecov.io](http://codecov.io/github/rafaqz/Dispersal.jl/coverage.svg?branch=master)](http://codecov.io/github/rafaqz/Dispersal.jl?branch=master)


This package builds on [Cellular.jl](https://github.com/rafaqz/Cellular.jl) to
provide a framework for cellular dispersal simulations.

To run a simulation:

```julia
using Cellular

# add a 0.0-1.0 scaled raster array that represents habitat suitability
suitability = your_2d_array

# define the model
model = Dispersal(layers=SuitabilityLayer(suitability))

# define the source array
source = zeros(Int8, size(suitability))

# seed it
source[24, 354] = 1

# run the simulation
sim!(source, model) 
```

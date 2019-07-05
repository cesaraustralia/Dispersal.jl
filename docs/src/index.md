# Dispersal.jl

```@meta
CurrentModule = Dispersal
```

```@docs
Dispersal
```

## Examples

This is the general pattern of running model in Dispersal.jl.

We define initial conditions. Then we run a dispersal
simulation that combines local and jump dispersal, for three timesteps.

```@example
using Cellular
using Dispersal

# Make an init array the same size as your suitability layer, and seed it
init = [0  0  0  0  0  0  0;
        0  0  0  0  0  0  0;
        0  0  0  0  0  0  0;
        0  0  0  0  0  0  0;
        0  0  1  0  0  0  0;
        0  0  0  0  0  0  0;
        0  0  0  0  0  0  0]

# add a 0.0-1.0 scaled raster array that represents habitat suitability
suitability = [0.5 0.0 0.3 0.0 0.0 0.3 0.0;
               0.0 0.2 0.8 0.0 0.9 0.6 0.0;
               0.0 0.5 1.0 1.0 1.0 0.0 0.0;
               0.0 0.0 1.0 1.0 0.7 0.0 0.0;
               0.0 0.0 0.6 1.0 0.0 0.0 0.0;
               0.0 0.1 0.6 0.0 0.0 0.0 0.0;
               0.0 0.0 0.0 0.0 0.0 0.0 0.7]


# Define a dispersal kernel function

f = d -> e^-d

# Define the neighborhood, using the dispersal kernel and a radius
hood = DispersalKernel(; f=f, radius=2, init=init)

# Define disersal models
localdisp = InwardsLocalDispersal(neighborhood=hood)
jumpdisp = JumpDispersal()
growth = SuitabilityExactLogisticGrowth(suitability)

# Set the output type
output = ArrayOutput(init)

# Run the simulation
sim!(output, Models(localdisp, jumpdisp, ), init, layers; time=1:3)

output.frames[3]
```

## Neighborhood Models

Models that consider the neighborhood of cells surrounding the current cell.
These disperse inwards to the current cell from the surrounding cell.

```@docs
AbstractInwardsDispersal
InwardsBinaryDispersal
InwardsPopulationDispersal
PoissonInwardsPopulationDispersal
```

## Partial Neighborhood Models

Neighborhood model that only operate on non-zero cells, dispersing outwards.

```@docs
AbstractOutwardsDispersal
OutwardsBinaryDispersal
OutwardsPopulationDispersal
```

### Neighborhoods

Kernel use neighborhoods, and extend Cellular.AbstractNeighborhood and `neighbors()` methods.

```@docs
AbstractDispersalKernel
DispersalKernel
DispersalKernel(; f=exponential, param=1.0, init=[], cellsize=1.0, radius=3)
```

## Cell Models

Models that simply transform the state of a single cell, ignoring the rest of the grid.


### Growth models

```@docs
AbstractGrowthModel
EulerExponentialGrowth
EulerLogisticGrowth
SuitabilityEulerExponentialGrowth
SuitabilityEulerLogisticGrowth
ExactExponentialGrowth
ExactLogisticGrowth
SuitabilityExactExponentialGrowth
SuitabilityExactLogisticGrowth
SuitabilityMask
```


### Mask layers

```@docs
Mask
```


### Allee effects

```@docs
AbstractAlleeExtinction
AlleeExtinction
```

## Partial Models

These models trigger rules that only partially update the grid.
The often operate only on cells that are currently occupied.

### Jump dispersal

```@docs
AbstractJumpDispersal
JumpDispersal
```

### Human driven dispersal

```@docs
AbstractHumanDispersal
HumanDispersal
```


## Layers

Layers provide overlay grids of raster data to models. They can be simple
matrices, or sequences for time series data.

Like models, than can be combined arbitrarily, in this case in a tuple. Methods
loop through all layers to return a scalar that is the product of their
outputs.

### Types

```@docs
AbstractSequence
Sequence
```

### Methods

```@docs
sequence_interpolate
cyclic
```

## Optimisation

Dispersal.jl provides optimisation tools for optimising the parameters of
arbitrary simulation rulesets. [`AbstractObjective`](@ref) can be extended
to add specific objection functions to transform simulation outputs.

```@docs
Parametriser
(p::Parametriser)(params)
AbstractObjective
SimpleObjective
RegionObjective
simpredictions
targets
```

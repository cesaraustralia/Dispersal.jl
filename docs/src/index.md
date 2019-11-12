# Dispersal.jl

```@meta
CurrentModule = Dispersal
```

```@docs
Dispersal
```

## Examples

This is the general pattern of running rules in Dispersal.jl.

We define initial conditions. Then we run a dispersal
simulation that combines local and jump dispersal, for three timesteps.

```@example
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


# Define the neighborhood, using the dispersal kernel and a radius
λ = 1.0
radius = 1
hood = DispersalKernel{radius}(;kernel=zeros(radius, radius), formulation=ExponentialKernel(λ))

# Define disersal rules
localdisp = InwardsBinaryDispersal(neighborhood=hood)

# Set the output type
output = ArrayOutput(init, 3)

# Run the simulation
sim!(output, Ruleset(localdisp); init=init, tspan=(1, 3))

output[3]
```

## Neighborhood Rules

Rules that consider the neighborhood of cells surrounding the current cell.
These disperse inwards to the current cell from the surrounding cell.

```@docs
AbstractInwardsDispersal
InwardsBinaryDispersal
InwardsPopulationDispersal
PoissonInwardsPopulationDispersal
```

## Partial Neighborhood Rules

Partial neighborhood rules that disperse outwards to the neighborhood 
when local populations exist.

```@docs
AbstractOutwardsDispersal
OutwardsBinaryDispersal
OutwardsPopulationDispersal
```

## Dispersal kernels 

Kernels extend `DynamicGrids.AbstractNeighborhood` and `neighbors()` methods.

```@docs
AbstractDispersalKernel
DispersalKernel
```

## Cell rules

Rules that simply transform the state of a single cell, ignoring the rest of the grid.


### Growth rules

```@docs
AbstractGrowthRule
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

## Partial Rules

These rules only partially update the grid. They often operate only on cells that
are currently occupied.

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

Layers provide overlay grids of raster data to rules. They can be simple
matrices, or sequences for time series data.

Like rules, than can be combined arbitrarily, in this case in a tuple. Methods
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

Dispersal.jl provides optimisation tools for automatically optimising 
the parameters of arbitrary rulesets given target data. [`AbstractObjective`](@ref) 
can be extended to add specific objection functions to transform simulation outputs.

```@docs
Parametriser
(p::Parametriser)(params)
AbstractObjective
SimpleObjective
RegionObjective
simpredictions
targets
```

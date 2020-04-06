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

# Define a simple init array the same size as your suitability layer, and seed it.

init = [0  0  0  0  0  0  0;
        0  0  0  0  0  0  0;
        0  0  0  0  0  0  0;
        0  0  0  0  0  0  0;
        0  0  1  0  0  0  0;
        0  0  0  0  0  0  0;
        0  0  0  0  0  0  0]

# Define an array that represents habitat suitability
suitability = [0.5 0.0 0.3 0.0 0.0 0.3 0.0;
               0.0 0.2 0.8 0.0 0.9 0.6 0.0;
               0.0 0.5 1.0 1.0 1.0 0.0 0.0;
               0.0 0.0 1.0 1.0 0.7 0.0 0.0;
               0.0 0.0 0.6 1.0 0.0 0.0 0.0;
               0.0 0.1 0.6 0.0 0.0 0.0 0.0;
               0.0 0.0 0.0 0.0 0.0 0.0 0.7]


# Define the neighborhood, using the dispersal kernel and a radius
hood = DispersalKernel{2}(; formulation=ExponentialKernel(1.0), 
                            distancemethod=CentroidToCentroid())

# Define disersal rules
localdisp = InwardsBinaryDispersal(neighborhood=hood)
jumpdisp = JumpDispersal()
growth = MaskGrowthMap(; layer=suitability)

# Set the output type
output = ArrayOutput(init, 3)

# Run the simulation
sim!(output, Ruleset(localdisp, jumpdisp, growth); init=init, tspan=(1, 3))

display(output[1])
display(output[2])
display(output[3])
```


## Neighborhood Rules

Rules that consider the neighborhood of cells surrounding the current cell.
These disperse population inwards to the current cell when populations exist 
in the surrounding cells.

```@docs
InwardsDispersal
InwardsBinaryDispersal
InwardsPopulationDispersal
PoissonInwardsPopulationDispersal
```

## Partial Neighborhood Rules

Partial neighborhood rules that disperse outwards to the neighborhood 
when a local population exists in the current cell. These methods
are harder to optimise and will generally have worse performance.

```@docs
OutwardsDispersal
OutwardsBinaryDispersal
OutwardsPopulationDispersal
```

## Dispersal kernels 

Kernels extend `DynamicGrids.Neighborhood`, and use `neighbors()` methods.

```@docs
DispersalKernel
```

### Distance methods

Dispersal kernels can be calculated in a number of ways, giving different 
properties and dispersal rates due to interactions between the cell size
and the dispersal distance.

```@docs
DistanceMethod
CentroidToArea
CentroidToCentroid
AreaToCentroid
AreaToArea
```

## Cell rules

Rules that simply transform the state of a single cell, ignoring the rest of the grid.


### Growth rules

```@docs
ExactExponentialGrowth
ExactLogisticGrowth
ExactExponentialGrowthMap
ExactLogisticGrowthMap
MaskGrowthMap
```


### Allee effects

```@docs
AlleeExtinction
```

## Partial Rules

These rules only partially update the grid. They often operate only on cells that
are currently occupied.

### Jump dispersal

```@docs
JumpDispersal
```

### Human driven dispersal

```@docs
HumanDispersal
populate!
populate
```


## Optimisation

Dispersal.jl provides optimisation tools for automatically optimising 
the parameters of arbitrary rulesets, given target data. [`Objective`](@ref) 
can be extended to add specific objection functions to transform simulation outputs.

```@docs
Parametriser
Objective
SimpleObjective
RegionObjective
ColorRegionFit
targets
predictions
```

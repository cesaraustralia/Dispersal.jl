# Dispersal.jl

```@meta
CurrentModule = Dispersal
```

```@docs
Dispersal
```

## Neighborhood Rules

Rules that consider the neighborhood of cells surrounding the current cell.
These disperse population inwards to the current cell when populations exist 
in the surrounding cells.

```@docs
InwardsDispersal
InwardsBinaryDispersal
InwardsPopulationDispersal
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
KernelFormulation
ExponentialKernel
```

### Distance methods

Dispersal kernels can be calculated in a number of ways, giving different 
properties and dispersal rates due to interactions between the cell size
and the dispersal distance.

```@docs
DistanceMethod
CentroidToCentroid
AreaToCentroid
AreaToArea
```

## Cell rules

Rules that simply transform the state of a single cell, ignoring the rest of the grid.


### Growth rules

```@docs
GrowthRule
ExactExponentialGrowth
ExactLogisticGrowth
GrowthMapRule
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
downsample
downsample!
```

# Layer

```@docs
layer
LayerCopy
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
RegionOutput
ColorRegionFit
Replicates
DistributedReplicates
SingleCoreReplicates
ThreadedReplicates
targets
predictions
```

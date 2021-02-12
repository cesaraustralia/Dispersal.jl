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
```

## Partial Neighborhood Rules

Partial neighborhood rules that disperse outwards to the neighborhood 
when a local population exists in the current cell. These methods
are harder to optimise and will generally have worse performance.

```@docs
OutwardsDispersal
```

## Dispersal kernels 

Kernels extend `DynamicGrids.Neighborhood`, and use `neighbors()` methods.

```@docs
DispersalKernel
KernelFormulation
ExponentialKernel
GeometricKernel
GaussianKernel
WeibullKernel
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
ExponentialGrowth
LogisticGrowth
ThresholdGrowth
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

### Mortality

```@docs
ExponentialMortality
LoglogisticMortality
```

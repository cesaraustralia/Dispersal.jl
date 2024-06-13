# Dispersal.jl

```@meta
CurrentModule = Dispersal
```

```@docs
Dispersal
```

## Growth rules

```@docs
GrowthRule
ExponentialGrowth
LogisticGrowth
ThresholdGrowth
```

## Mortality

```@docs
Mortality
ExponentialMortality
LoglogisticMortality
```

## Allee effects

```@docs
AlleeExtinction
```

## Local dispersal rules

```@docs
InwardsDispersal
OutwardsDispersal
```

### Dispersal kernels 

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

## Jump dispersal

```@docs
JumpDispersal
```

## Human driven dispersal

```@docs
HumanDispersal
populate!
populate
```

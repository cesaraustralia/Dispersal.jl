module Dispersal
# Use the README as the module docs
@doc read(joinpath(dirname(@__DIR__), "README.md"), String) Dispersal

using ConstructionBase,
      Colors,
      Dates,
      DimensionalData,
      Distributed,
      Distributions,
      DocStringExtensions,
      Flatten,
      LinearAlgebra,
      LossFunctions,
      ModelParameters,
      PoissonRandom,
      Reexport,
      Setfield,
      Statistics

@reexport using DynamicGrids
@reexport using ModelParameters

using LossFunctions: ZeroOneLoss

import Base: getindex, setindex!, lastindex, size, length, push!

import DynamicGrids: applyrule, applyrule!, precalcrule,
       neighbors, sumneighbors, neighborhood, setneighbor!, mapsetneighbor!,
       radius, gridsize, mask, overflow, cellsize, ruleset, inbounds, extent,
       currenttime, currenttimestep, timestep, tspan,
       buffer, SimData, WritableGridData, aux, unwrap, kernel

import ConstructionBase: constructorof


export AbstractDispersalKernel, DispersalKernel

export KernelFormulation, ExponentialKernel, GeometricKernel, GaussianKernel

export DistanceMethod, CentroidToCentroid, CentroidToArea, AreaToArea, AreaToCentroid

export AbstractInwardsDispersal, InwardsBinaryDispersal, InwardsPopulationDispersal,
       SwitchedInwardsPopulationDispersal

export AbstractOutwardsDispersal, OutwardsBinaryDispersal, OutwardsPopulationDispersal

export AlleeExtinction, JumpDispersal, HumanDispersal

export BatchGroups, HeirarchicalGroups

export ExponentialGrowth, LogisticGrowth

export DiscreteGrowth, DiscreteGrowthMap, DiscreteGrowth2

export MaskGrowthMap, ExactExponentialGrowthMap, ExactLogisticGrowthMap, ExactLogisticGrowthMap2, ExactLogisticGrowthMap3

export SurvLogLogistic, SurvLogLogisticMap, SurvLogLogisticMap2

export GrowthSurvLogLogisticMap3, SurvLogLogisticMapTuple, GrowthMapTuple

export MaskSurvMap

export SelectionGradientSurv, SelectionGradientSurvMap

export SelectionGradient1locusSurv, SelectionGradient1locusSurvMap, LandeVariable, SelectionGradientMapTuple

export DeltaAlleleFrequencySurv, DeltaAlleleFrequencySurvMap, DeltaAlleleFrequencySurvMap_noFastmath

export MatingPopulation, MatingDispersal

export MaskGrowthMap, ExponentialGrowthMap, LogisticGrowthMap

export AuxCopy

export Parametriser, AbstractObjective, SimpleObjective, RegionObjective, RegionOutput,
       ColorRegionFit, Accuracy

export ThreadedReplicates, DistributedReplicates, SingleCoreReplicates

export targets, predictions


# Documentation templates
@template TYPES =
    """
    $(TYPEDEF)
    $(DOCSTRING)
    """

include("downsampling.jl")
include("utils.jl") 
include("rules/kernel/common.jl")
include("rules/kernel/inwards.jl")
include("rules/kernel/outwards.jl")
include("rules/growth.jl")
include("rules/human.jl")
include("rules/jump.jl")
include("rules/allee.jl")
include("rules/survival.jl")
include("rules/allele_frequency.jl")
include("rules/selection_gradient.jl")
include("rules/mating_population.jl")
include("rules/discrete_growth.jl")
include("rules/auxcopy.jl")
include("optimisation/optimisation.jl")
include("optimisation/objectives.jl")
include("optimisation/output.jl")


end # module

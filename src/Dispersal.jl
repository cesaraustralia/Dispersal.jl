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
      FieldDocTables,
      Flatten,
      LinearAlgebra,
      LossFunctions,
      PoissonRandom,
      Reexport,
      Setfield,
      Statistics

@reexport using DynamicGrids
@reexport using ModelParameters

using LossFunctions: ZeroOneLoss
using DimensionalData: Time

import Base: getindex, setindex!, lastindex, size, length, push!

import DynamicGrids: applyrule, applyrule!, precalcrule,
       neighbors, sumneighbors, neighborhood, setneighbor!, mapsetneighbor!,
       radius, gridsize, mask, overflow, cellsize, ruleset, inbounds, extent,
       currenttime, currenttimestep, timestep, tspan,
       buffer, SimData, WritableGridData, aux, unwrap

import ConstructionBase: constructorof


export AbstractDispersalKernel, DispersalKernel

export KernelFormulation, ExponentialKernel

export DistanceMethod, CentroidToCentroid, CentroidToArea, AreaToArea, AreaToCentroid

export AbstractInwardsDispersal, InwardsBinaryDispersal, InwardsPopulationDispersal,
       SwitchedInwardsPopulationDispersal

export AbstractOutwardsDispersal, OutwardsBinaryDispersal, OutwardsPopulationDispersal


export AlleeExtinction, JumpDispersal, HumanDispersal

export BatchGroups, HeirarchicalGroups

export ExactExponentialGrowth, ExactLogisticGrowth

export MaskGrowthMap, ExactExponentialGrowthMap, ExactLogisticGrowthMap

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

const FIELDDOCTABLE = FieldDocTable((;))


include("downsampling.jl")
include("utils.jl") 
include("rules/kernel/common.jl")
include("rules/kernel/inwards.jl")
include("rules/kernel/outwards.jl")
include("rules/growth.jl")
include("rules/human.jl")
include("rules/jump.jl")
include("rules/allee.jl")
include("rules/auxcopy.jl")
include("optimisation/optimisation.jl")
include("optimisation/objectives.jl")
include("optimisation/output.jl")


end # module

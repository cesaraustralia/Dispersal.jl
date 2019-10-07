module Dispersal
# Use the README as the module docs
@doc read(joinpath(dirname(@__DIR__), "README.md"), String) Dispersal

using Colors,
      DimensionalData,
      Distributed,
      DocStringExtensions,
      FieldDefaults,
      FieldDocTables,
      FieldMetadata,
      Flatten,
      LinearAlgebra,
      LossFunctions,
      Mixers,
      PoissonRandom,
      Reexport,
      Setfield,
      Statistics

@reexport using DynamicGrids

using LossFunctions: ZeroOneLoss
using DimensionalData: Time
using Interpolations: WeightedArbIndex

import Base: getindex, setindex!, lastindex, size, length, push!


import DynamicGrids: applyrule, applyrule!, neighbors, radius,
       currenttime, framesize, mask, overflow, timestep, cellsize, ruleset


import Flatten: constructor_of

import FieldMetadata: @description, @limits, @flattenable,
                      default, description, limits, flattenable



export AbstractDispersal

export AbstractInwardsDispersal, InwardsBinaryDispersal, InwardsPopulationDispersal,
       PoissonInwardsPopulationDispersal

export AbstractOutwardsDispersal, OutwardsBinaryDispersal, OutwardsPopulationDispersal

export AbstractDispersalKernel, DispersalKernel, AbstractKernelFormulation,
       ExponentialKernel

export AbstractDistanceMethod, CentroidToCentroid, CentroidToArea, AreaToArea, AreaToCentroid

export AbstractJumpDispersal, JumpDispersal

export AbstractHumanDispersal, HumanDispersal #, populate, populate!

export ExactExponentialGrowth, ExactLogisticGrowth

export MaskGrowthMap, ExactExponentialGrowthMap, ExactLogisticGrowthMap

export AlleeExtinction, AlleeExtinction

export Parametriser, AbstractObjective, SimpleObjective, RegionObjective, RegionOutput, 
       ColorRegionFit, Accuracy

export ThreadedReplicates, DistributedReplicates, SingleCoreReplicates

export targets, predictions


const FIELDDOCTABLE = FieldDocTable((:Description, :Default, :Limits),
                                    (description, default, limits);
                                    truncation=(100,40,100))

# Documentation templates
@template TYPES =
    """
    $(TYPEDEF)
    $(DOCSTRING)
    """

include("downsampling.jl")
include("layers.jl")
include("rules/mixins.jl")
include("rules/kernel/common.jl")
include("rules/kernel/inwards.jl")
include("rules/kernel/outwards.jl")
include("rules/growth.jl")
include("rules/human.jl")
include("rules/jump.jl")
include("rules/allee.jl")
include("optimisation/optimisation.jl")
include("optimisation/objectives.jl")
include("optimisation/frame_processing.jl")


end # module

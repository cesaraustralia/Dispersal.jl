module Dispersal
# Use the README as the module docs
@doc read(joinpath(dirname(@__DIR__), "README.md"), String) Dispersal

using ConstructionBase,
      Dates,
      DimensionalData,
      Distributions,
      DocStringExtensions,
      LinearAlgebra,
      ModelParameters,
      Reexport,
      Setfield,
      StaticArrays

@reexport using DynamicGrids, 
                ModelParameters

using ModelParameters: params

import DynamicGrids: applyrule, applyrule!, precalcrule, kernel

export KernelFormulation, ExponentialKernel

export DistanceMethod, CentroidToCentroid, AreaToArea, AreaToCentroid #, CentroidToArea, 

export DispersalKernel

export InwardsDispersal, OutwardsDispersal

export AlleeExtinction, JumpDispersal

export HumanDispersal, BatchGroups, HeirarchicalGroups

export ThresholdGrowth, ExponentialGrowth, LogisticGrowth

export LogisticSurvival

export Pulsed_Exposure, Degradation_Exposure, Threshold_Exposure

export Lande_Resistance

# Documentation templates
@template TYPES =
    """
    $(TYPEDEF)
    $(DOCSTRING)
    """

const DG = DynamicGrids

include("downsampling.jl")
include("utils.jl") 
include("kernel/common.jl")
include("kernel/inwards.jl")
include("kernel/outwards.jl")
include("growth.jl")
include("human.jl")
include("jump.jl")
include("allee.jl")
include("exposure.jl")
include("survival.jl")
include("resistance.jl")

end # module

module Dispersal
# Use the README as the module docs
@doc read(joinpath(dirname(@__DIR__), "README.md"), String) Dispersal

using Dates,
      Distributions,
      DocStringExtensions,
      LinearAlgebra,
      Reexport

@reexport using DynamicGrids 

using DynamicGrids.ConstructionBase
using DynamicGrids.DimensionalData
using DynamicGrids.ModelParameters
using DynamicGrids.Setfield
using DynamicGrids.Stencils
using DynamicGrids.StaticArrays
using Stencils

using DynamicGrids.ModelParameters: params

import DynamicGrids: applyrule, applyrule!, modifyrule, kernel

export KernelFormulation, ExponentialKernel, GeometricKernel, GaussianKernel, WeibullKernel

export DistanceMethod, CentroidToCentroid, AreaToArea, AreaToCentroid #, CentroidToArea, 

export DispersalKernel

export InwardsDispersal, OutwardsDispersal

export AlleeExtinction, JumpDispersal

export HumanDispersal, BatchGroups, HeirarchicalGroups

export ExponentialGrowth, LogisticGrowth, ThresholdGrowth

export LoglogisticMortality, ExponentialMortality

# Documentation templates
@template TYPES =
    """
    $(TYPEDEF)
    $(DOCSTRING)
    """

const DG = DynamicGrids

include("downsampling.jl")
include("kernel/formulations.jl")
include("kernel/distancemethods.jl")
include("kernel/kernel.jl")
include("kernel/inwards.jl")
include("kernel/outwards.jl")
include("growth.jl")
include("human.jl")
include("jump.jl")
include("allee.jl")
include("mortality.jl")

end # module

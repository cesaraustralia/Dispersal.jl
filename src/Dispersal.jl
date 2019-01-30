"""
This package extends [Cellular.jl](https://github.com/rafaqz/Cellular.jl)

[Dispersal.jl](https://github.com/rafaqz/Dispersal.jl) provides a range of
dispersal modules that can be combined to build grid-based organism dispersal simulations.

The framework is highly extensible. A model may start with the defaults and formulations provided,
but incrementally customise them for a particular use-case.

Additionally, components provided by Cellular.jl or other
packages may be incorporated into a simulaiton.
"""
module Dispersal

using Cellular,
      DocStringExtensions,
      LinearAlgebra,
      Parameters,
      Mixers,
      Flatten,
      Requires,
      FieldMetadata,
      Distributions,
      Distributed,
      Statistics,
      Colors


import Base: getindex, setindex!, lastindex, size, length, push!

import Cellular: rule, rule!, neighbors, inbounds, radius, temp_neighborhood

import FieldMetadata: @description, @limits, @flattenable, description, limits, flattenable



export AbstractDispersal

export AbstractInwardsDispersal, InwardsBinaryDispersal, InwardsPopulationDispersal

export AbstractOutwardsDispersal, OutwardsBinaryDispersal, OutwardsPopulationDispersal

export AbstractJumpDispersal, JumpDispersal

export AbstractHumanDispersal, HumanDispersal, populate, populate!

export EulerExponentialGrowth, EulerLogisticGrowth, ExactExponentialGrowth, ExactLogisticGrowth

export SuitabilityMask, SuitabilityEulerExponentialGrowth, SuitabilityEulerLogisticGrowth,
       SuitabilityExactExponentialGrowth, SuitabilityExactLogisticGrowth

export AlleeExtinction, AlleeExtinction

export Mask

export AbstractDispersalKernel, DispersalKernel, exponential 
 
export Sequence

export RegionParametriser, SumOutput, ColorRegionFit


# Documentation templates
@template TYPES =
    """
    $(TYPEDEF)
    $(FIELDS)
    """


include("downsampling.jl")
include("layers.jl")
include("models/mixins.jl")
include("models/kernel/common.jl")
include("models/kernel/inwards.jl")
include("models/kernel/outwards.jl")
include("models/growth.jl")
include("models/mask.jl")
include("models/human.jl")
include("models/jump.jl")
include("models/allee.jl")
include("optimisation/region.jl")
include("optimisation/cell.jl")


end # module

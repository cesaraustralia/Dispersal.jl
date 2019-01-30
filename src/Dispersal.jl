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
      Statistics

import Cellular: rule, rule!, neighbors, inbounds, radius, temp_neighborhood
import Base: getindex, setindex!, lastindex, size, length, push!
import Flatten: flattenable
import FieldMetadata: @description, @limits, description, limits

# Documentation templates
@template TYPES =
    """
    $(TYPEDEF)
    $(FIELDS)
    """

include("types.jl")
include("downsampling.jl")
include("layers.jl")
include("kernel/common.jl")
include("kernel/inwards.jl")
include("kernel/outwards.jl")
include("growth.jl")
include("mask.jl")
include("human.jl")
include("jump.jl")
include("allee.jl")

export AbstractDispersal,
       AbstractInwardsDispersal,
       InwardsBinaryDispersal,
       InwardsPopulationDispersal,
       AbstractOutwardsDispersal,
       OutwardsBinaryDispersal,
       OutwardsPopulationDispersal,
       HudginsDispersal,
       AbstractJumpDispersal,
       JumpDispersal,
       AbstractHumanDispersal,
       HumanDispersal,
       EulerExponentialGrowth,
       EulerLogisticGrowth,
       ExactExponentialGrowth,
       ExactLogisticGrowth,
       AlleeExtinction,
       Mask,
       SuitabilityMask,
       SuitabilityEulerExponentialGrowth,
       SuitabilityEulerLogisticGrowth,
       SuitabilityExactExponentialGrowth,
       SuitabilityExactLogisticGrowth,
       AbstractDispersalKernel,
       DispersalKernel,
       DispersalGrid,
       Layer,
       Sequence,
       exponential,
       populate,
       populate!


end # module

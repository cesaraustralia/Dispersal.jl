"""
This package extends [Cellular.jl](https://github.com/rafaqz/Cellular.jl)

[Dispersal.jl](https://github.com/rafaqz/Dispersal.jl) provides a range of
dispersal modules that can be combined to build grid-based organism dispersal simulations.

The design provides a solid framework while allowing customisation of any 
aspect. A model may start with the defaults and formulations provided, 
but incrementally customise them for a particular use-case, to any level of detail. 

Additionally, modules, outputs, neighborhoods provided by Cellular.jl or other 
packages that extend it may be incorporated into a simulaiton.
"""
module Dispersal

using Cellular, 
      DocStringExtensions, 
      Parameters, 
      Mixers, 
      Flatten, 
      Requires,
      FieldMetadata

import Cellular: rule, rule!, neighbors, inbounds
import Base: getindex, setindex!, endof, size, length, push!
import Flatten: flattenable
import FieldMetadata: @description, @limits, description, limits

# Documentation templates
@template TYPES =
    """
    $(TYPEDEF)
    $(FIELDS)
    """

include("types.jl")
include("utils.jl")
include("kernels.jl")
include("neighborhoods.jl")
include("rules.jl")
include("layers.jl")
include("growth.jl")
include("human.jl")

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
       ExponentialGrowth,
       LogisticGrowth,
       SuitabilityMask, 
       SuitabilityExponentialGrowth,
       SuitabilityLogisticGrowth,
       AbstractDispersalNeighborhood, 
       DispersalNeighborhood,
       DispersalGrid,
       Sequencwe,
       Layer,
       SuitabilityLayer, 
       SuitabilitySequence,
       HumanLayer, 
       exponential

function __init__()
    @require CuArrays="3a865a2d-5b23-5a0f-bc46-62713ec82fae" begin
        include("cuda.jl")
        include("hudgins.jl")
    end
end



end # module

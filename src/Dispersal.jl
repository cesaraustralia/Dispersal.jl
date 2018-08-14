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

using Cellular, DocStringExtensions, Parameters, Mixers

import Cellular: rule, neighbors, inbounds
import Base: getindex, setindex!, endof, size, length, push!

# Documentation templates
@template TYPES =
    """
    $(TYPEDEF)
    $(FIELDS)
    """

include("types.jl")
include("neighborhoods.jl")
include("rules.jl")
include("layers.jl")

export AbstractDispersal,
       AbstractInwardsDispersal, 
       InwardsLocalDispersal,
       AbstractOutwardsDispersal, 
       OutwardsLocalDispersal,
       HudginsDispersal,
       AbstractJumpDispersal, 
       JumpDispersal,
       AbstractHumanDispersal, 
       HumanDispersal,
       AbstractDispersalNeighborhood, 
       DispersalNeighborhood,
       DispersalGrid,
       AbstractLayers,
       Layers,
       AbstractLayer, 
       AbstractSuitabilityLayer, 
       SuitabilityLayer, 
       Layer,
       AbstractHumanLayer,
       HumanLayer,
       AbstractSuitabilitySequence, 
       SuitabilitySequence,
       exponential

end # module

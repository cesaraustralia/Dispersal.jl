"""
This package extends [Cellular.jl](https://github.com/rafaqz/Cellular.jl)

[Dispersal.jl](https://github.com/rafaqz/Dispersal.jl) provides a range of
dispersal modules that can be combined to build grid-based organism dispersal simulations.

The design provides a solid framewhork while allowing customisation of any 
aspect of a simulation. A model may start with the defaults and formulations provided, 
but incrementally customise them for a particular use-case, to any level of detail. 

Additionally, modules, outputs, neighborhoods provided by Cellular.jl or other 
packages that extend it may all be incorporated into a simulaiton.
"""
module Dispersal

using Cellular, DocStringExtensions, Parameters

import Cellular: rule, neighbors

include("disperse.jl")
include("layers.jl")

export AbstractDispersal,
       AbstractLocalDispersal, 
       LocalDispersal,
       AbstractJumpDispersal, 
       JumpDispersal,
       AbstractHumanDispersal, 
       HumanDispersal,
       AbstractDispersalNeighborhood, 
       DispersalNeighborhood,
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

module Dispersal

using Cellular
using Parameters
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
       HumanLayer, 
       AbstractSuitabilityLayer, 
       SuitabilityLayer,
       AbstractSuitabilitySequence, 
       SuitabilitySequence,
       exponential

end # module

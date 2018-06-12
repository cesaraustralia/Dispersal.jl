module Dispersal

using Cellular
using Parameters

include("dispersal.jl")
include("layers.jl")

export AbstractDispersal,
       StratifiedDispersal,
       AbstractShortDispersal, 
       ShortDispersal,
       AbstractLongDispersal, 
       JumpDispersal,
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

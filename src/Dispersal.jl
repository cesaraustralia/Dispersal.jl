"""
[Dispersal.jl](https://github.com/rafaqz/Dispersal.jl) extends
[Cellular.jl](https://github.com/rafaqz/DynamicGrids.jl) for grid-based organism dispersal simulations.

Dispersal.jl provides a range of simulation rules than can be combined to develop complex organism
dispersal models. Growth rates, dispersal kernels, allee effects, and randomised jump
and human assisted dispersal rules are provided.

These rules can be chained arbitrarily, and custom rules combined with the provided set.

DynamicGrids.jl includes REPLOutput for live display in the REPL, while
[CellularAutomataGtk](https://github.com/rafaqz/CellularAutomataGtk.jl) provides GtkOutput for a
simple graphical viewer, and [CellularAutomataWeb](https://github.com/rafaqz/CellularAutomataWeb.jl)
provides BlinkOuput and MuxServer for desktop and online web applications.
These web apps automatically provide realtime slider controls for model parameters,
including for custom rules, customisable using [FieldMetadata.jl](https://github.com/rafaqz/FieldMetadata.jl).

[GrowthRates.jl](https://github.com/rafaqz/GrowthRates.jl) can efficiently generate
the layers required for suitability growth rules based on temperature response and stress factors.
"""
module Dispersal

using Colors,
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
      Statistics

@reexport using DynamicGrids

using LossFunctions: ZeroOneLoss

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

export EulerExponentialGrowth, EulerLogisticGrowth, ExactExponentialGrowth,
       ExactLogisticGrowth

export SuitabilityMask, SuitabilityEulerExponentialGrowth, SuitabilityEulerLogisticGrowth,
       SuitabilityExactExponentialGrowth, SuitabilityExactLogisticGrowth

export AlleeExtinction, AlleeExtinction

export Sequence

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

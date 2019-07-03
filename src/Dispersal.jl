"""
[Dispersal.jl](https://github.com/rafaqz/Dispersal.jl) extends
[Cellular.jl](https://github.com/rafaqz/Cellular.jl) to provides a range of
dispersal modules that can be combined to build grid-based organism dispersal simulations.

Dispersal.jl aims to provide a set of rules than can be combined to develope complex organidm dispersal models.
Growth rates, dispersal kernels, allee effects, and randomised jump and human assisted dispersal modes are provided.

Rules can be chained arbitrarily, and cutom rules mixed with the provided set.
The framework is flexible and extensible, and basically anything can be customised,
from outputs, rules and optimisation methods.

Outputs from Cellular.jl include REPLOutput for live display in the REPL, Gtk output for a
simple graphical window, and BlinkOuput and MuxServer for automatic desktop and served web applications.
Web apps automatically provide realtime slider controls for model parameters, including custom rules.

[GrowthRates.jl](https://github.com/rafaqz/GrowthRates.jl) can efficiently generate
the layers required for suitability growth rules based on temperature response and stress factors.
"""
module Dispersal

using CellularAutomataBase,
      Colors,
      Distributed,
      Distributions,
      DocStringExtensions,
      FieldDefaults,
      FieldDocTables,
      FieldMetadata,
      Flatten,
      LinearAlgebra,
      LossFunctions,
      Mixers,
      Statistics

using LossFunctions: ZeroOneLoss

import Base: getindex, setindex!, lastindex, size, length, push!

import CellularAutomataBase: applyrule, applyrule!, neighbors, timestep, radius, buffer

import FieldMetadata: @description, @limits, @flattenable, 
                      default, description, limits, flattenable


export AbstractDispersal

export AbstractInwardsDispersal, InwardsBinaryDispersal, InwardsPopulationDispersal,
       PoissonInwardsPopulationDispersal

export AbstractOutwardsDispersal, OutwardsBinaryDispersal, OutwardsPopulationDispersal

export AbstractDispersalKernel, DispersalKernel, AbstractKernelFormulation, 
       ExponentialKernel 

export AbstractJumpDispersal, JumpDispersal

export AbstractHumanDispersal, HumanDispersal, populate, populate!

export EulerExponentialGrowth, EulerLogisticGrowth, ExactExponentialGrowth, 
       ExactLogisticGrowth

export SuitabilityMask, SuitabilityEulerExponentialGrowth, SuitabilityEulerLogisticGrowth,
       SuitabilityExactExponentialGrowth, SuitabilityExactLogisticGrowth

export AlleeExtinction, AlleeExtinction

export Mask

export Sequence

export Parametriser, AbstractObjective, RegionObjective, ColorRegionFit



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
include("rules/mask.jl")
include("rules/human.jl")
include("rules/jump.jl")
include("rules/allee.jl")
include("optimisation.jl")


end # module

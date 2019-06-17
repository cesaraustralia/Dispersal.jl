"""
[Dispersal.jl](https://github.com/rafaqz/Dispersal.jl) extends
[Cellular.jl](https://github.com/rafaqz/Cellular.jl) to provides a range of
dispersal modules that can be combined to build grid-based organism dispersal simulations.

Dispersal.jl aims to provide a set of models than can be combined to develope complex organidm dispersal models.
Growth rates, dispersal kernels, allee effects, and randomised jump and human assisted dispersal modes are provided.

Models can be chained arbitrarily, and cutom models mixed with the provided set.
The framework is flexible and extensible, and basically anything can be customised,
from outputs, models and optimisation methods.

Outputs from Cellular.jl include REPLOutput for live display in the REPL, Gtk output for a
simple graphical window, and BlinkOuput and MuxServer for automatic desktop and served web applications.
Web apps automatically provide realtime slider controls for model parameters, including custom models.

[GrowthRates.jl](https://github.com/rafaqz/GrowthRates.jl) can efficiently generate
the layers required for suitability growth models based on temperature response and stress factors.
"""
module Dispersal

using Cellular,
      DocStringExtensions,
      FielddocTables,
      LinearAlgebra,
      FieldDefaults,
      Mixers,
      Flatten,
      Requires,
      FieldMetadata,
      Distributions,
      Distributed,
      Statistics,
      Colors

using LossFunctions: ZeroOneLoss

import Base: getindex, setindex!, lastindex, size, length, push!

import Cellular: rule, rule!, neighbors, inbounds, radius, temp_neighborhood

import FieldMetadata: @description, @limits, @flattenable, default, description, limits, flattenable


export AbstractDispersal

export AbstractInwardsDispersal, InwardsBinaryDispersal, InwardsPopulationDispersal,
       PoissonInwardsPopulationDispersal

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

export Parametriser, RegionObjective, SumOutput, ColorRegionFit

const FIELDDOCTABLE = FielddocTable((:Description, :Default, :Limits),
                                    (description, default, limits);
                                    truncation=(100,40,100))

# Documentation templates. Currently broken in macros, so basically everywhere in this package...
@template TYPES =
    """
    $(TYPEDEF)
    $(DOCSTRING)
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

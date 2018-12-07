"""
Layers wrap an array, normally built from some kind of raster.

The wrapper defines the purpose of the layer and allows specialised
method dispatch to utilise them.
"""
abstract type AbstractLayerMatrix{T,S<:AbstractMatrix{T}} <: AbstractMatrix{T} end
abstract type AbstractLayer{T,S} <: AbstractLayerMatrix{T,S} end

(::Type{F})(data::D) where {F<:AbstractLayer,D} = F{eltype(D),D}(data)

length(l::AbstractLayer) = length(l.data)
size(l::AbstractLayer) = size(l.data)
firstindex(l::AbstractLayer) = firstindex(l.data)
lastindex(l::AbstractLayer) = lastindex(l.data)
Base.@propagate_inbounds getindex(l::AbstractLayer, i...) = getindex(l.data, i...)
Base.@propagate_inbounds setindex!(l::AbstractLayer, x, i...) = setindex!(l.data, x, i...)
push!(l::AbstractLayer, x) = push!(l.data, x)

@premix struct Layer{T,S}
    "Any 2-dimensional AbstractArray matching the coordinates of the init array"
    data::S
end


abstract type AbstractHumanLayer{T,S} <: AbstractLayer{T,S} end

"A wrapper for arrays that provide human dispersal scalars for the grid."
@Layer struct HumanLayer{} <: AbstractHumanLayer{T,S} end


"Abstract type for layers that provide suitability scalars."
abstract type AbstractSuitabilityLayer{T,S} <: AbstractLayer{T,S} end

"A wrapper for arrays that provide suitability scalars for the grid."
@Layer struct SuitabilityLayer{} <: AbstractSuitabilityLayer{T,S} end


"""
An abstract type for for sequences of suitability layers.
"""
abstract type AbstractSequence{T} <: AbstractVector{T} end

length(s::AbstractSequence) = length(s.data)
size(s::AbstractSequence) = size(s.data, 1)
firstindex(s::AbstractSequence) = firstindex(s.data)
lastindex(s::AbstractSequence) = lastindex(s.data)
Base.@propagate_inbounds getindex(s::AbstractSequence{T}, i...) where T = 
    getindex(s.data, i...)


@premix struct Sequence{T,D,TS}
    "Any 2-dimensional AbstractArray matching the coordinates of the init array"
    data::D
    "The timespan of each layer in the sequence."
    timespan::TS
end

(::Type{F})(data::D, timespan::TS) where {F<:AbstractSequence,D,TS} = 
    F{eltype(D),D,TS}(data, timespan)

abstract type AbstractSuitabilitySequence{T} <: AbstractSequence{T} end

"""
A wrapper for a sequence of arrays that provide suitability scalar
for grid points that cycle over a time span.
"""
@Sequence struct SuitabilitySequence{} <: AbstractSuitabilitySequence{T} end


"""
    suitability(layers, index, t::Number)
Returns a scalar representing cell suitability from one or multiple layers.

For multiple layers the product of each individual scalar is returned.

Layers of type other than AbstractSuitabilityLayer return 1.0.

### Arguments
- `layers` : a single layer or tuple of layers of any type
- `index` : a tuple containing 2 grid coordinates
- `t::Number` : current timestep for interploating layere sequences.
"""
function suitability end
suitability(layers::Tuple, index, t::Number) = 
    suitability(layers[1], index, t) * suitability(Base.tail(layers), index, t::Number)
suitability(layers::Tuple{}, index, t::Number) = 1
suitability(layer::AbstractSuitabilityLayer, index, t::Number) = 
    @inbounds return layer[index...]
suitability(sequence::AbstractSuitabilitySequence, index, t::Number) =
    sequence_interpolate(sequence, index, t)
suitability(layers, index, t::Number) = 1

"""
    human_impact(layers, index, t)
Returns a scalar indicating human impact from one or multiple layers.

For multiple layers the product of each individual scalar is returned.

Layers of type other than AbstractHumanLayer return 1.0.

### Arguments
- `layers` : a single layer or tuple of layers of any type
- `index` : a tuple containing 2 grid coordinates
- `t::Number` : current timestep for interploating layere sequences.
"""
function human_impact end
human_impact(layers::Tuple, index, t::Number) = 
    human_impact(layers[1], index, t) * human_impact(Base.tail(layers), index, t::Number)
human_impact(layers::Tuple{}, index, t::Number) = 1
human_impact(layer::AbstractHumanLayer, index, t::Number) = @inbounds return layer[index...]
human_impact(layers, index, t::Number) = 1

"""
    sequence_interpolate(layer, index, t)
Interpolates between layers in a sequence.
"""
@inline sequence_interpolate(seq, index, t::Number) = begin
    # Time position is centered in the current frame, not at the start.
    tf = (t + seq.timespan / 2) / seq.timespan
    # Linear interpolation between layers in the sequence.
    tr = unsafe_trunc(Int64, tf)
    frac = tf - tr
    pos1 = cyclic(tr, length(seq))
    pos2 = cyclic(tr + 1, length(seq))
    @inbounds seq[pos1][index...] * (oneunit(eltype(seq[1])) - frac) + seq[pos2][index...] * frac
end

"Cycles a time position through a particular timespan length"
cyclic(t::Int, len::Int) =
    if t > len
        rem(t + len - 1, len) + 1
    elseif t <= 0
        t + (t ÷ len -1) * -len
    else
        t
    end

"""
Layers wrap an array, normally built from some kind of raster.

The wrapper defines the purpose of the layer and allows specialised
method dispatch to utilise them.
"""
abstract type AbstractLayerMatrix{T, S <: AbstractMatrix{T}} <: AbstractMatrix{T} end
abstract type AbstractLayer{T,S} <: AbstractLayerMatrix{T,S} end

(::Type{F})(data::D) where {F<:AbstractLayer,D} = F{eltype(D),D}(data)

length(l::AbstractLayer) = length(l.data)
size(l::AbstractLayer) = size(l.data)
endof(l::AbstractLayer) = endof(l.data)
getindex(l::AbstractLayer, i...) = getindex(l.data, i...)
setindex!(l::AbstractLayer, x, i...) = setindex!(l.data, x, i...)
push!(l::AbstractLayer, x) = push!(l.data, x)

abstract type AbstractHumanLayer{T,S} <: AbstractLayer{T,S} end

@premix struct Data{T,S<:AbstractMatrix}
    "Any 2-dimensional AbstractArray matching the coordinates of the init array"
    data::S
end

"A wrapper for arrays that provide human dispersal scalars for the grid."
@Data struct HumanLayer{} <: AbstractHumanLayer{T,S} end

"Abstract type for layers that provide suitability scalars."
abstract type AbstractSuitabilityLayer{T,S} <: AbstractLayer{T,S} end

"A wrapper for arrays that provide suitability scalars for the grid."
@Data struct SuitabilityLayer{} <: AbstractSuitabilityLayer{T,S} end

"""
An abstract type for for sequences of suitability layers.
"""
abstract type AbstractSequence{T} <: AbstractVector{T} end

length(s::AbstractSequence) = length(s.data)
size(s::AbstractSequence) = size(s.data)
endof(s::AbstractSequence) = endof(s.data)
getindex(s::AbstractSequence{T}, i) where T = getindex(s.data::T, i)::eltype(T)
setindex!(s::AbstractSequence, x, i) = setindex!(s.data, x, i)
push!(s::AbstractSequence, x) = push!(s.data, x)

abstract type AbstractSuitabilitySequence{T} <: AbstractSequence{T} end

"""
A wrapper for a sequence of arrays that provide suitability scalar
for grid points that cycle over a time span.
"""
struct SuitabilitySequence{T,S} <: AbstractSuitabilitySequence{T}
    "Any 2-dimensional AbstractArray matching the coordinates of the init array"
    data::T
    "The timespan of each layer in the sequence."
    timespan::S
end

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
suitability(layers, index, t::Number) = 1.0
suitability(layers::Tuple, index, t::Number) =
    mapreduce(l -> suitability(l, index, t), *, 1.0, layers)
suitability(layer::AbstractSuitabilityLayer, index, t::Number) = layer[index...]
suitability(sequence::AbstractSuitabilitySequence, index, t::Number) =
    sequence_interpolate(sequence, index, t)


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
human_impact(layers::Tuple, index, t) =
    mapreduce(l -> human_impact(l, index, t), *, 1.0, layers)
human_impact(layer::AbstractHumanLayer, index, t) = layer[round(Int,t)][index...]
human_impact(layer, index, t) = 1.0

"""
    sequence_interpolate(layer, index, t)
Interpolates between layers in a sequence.
"""
sequence_interpolate(sequence, index, t::Number) = begin
    # Time position is centered in the current frame, not at the start.
    tf = (t + sequence.timespan * 0.5) / sequence.timespan
    # Linear interpolation between layers in the sequence.
    tr = floor(Int, tf)
    frac = convert(eltype(sequence[1]), tf - tr)
    pos1 = cyclic(tr, length(sequence))
    pos2 = cyclic(tr + 1, length(sequence))
    sequence[pos1][index[1], index[2]] * (oneunit(eltype(sequence[1])) - frac) + sequence.data[pos2][index...] * frac
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

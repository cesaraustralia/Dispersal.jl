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
show(s::AbstractSequence) = show(s.data)

struct Sequence{T,TS} <: AbstractSequence{T}
    "Any 2-dimensional AbstractArray matching the coordinates of the init array"
    data::T
    "The length in time of each layer in the sequence."
    timestep::TS
end


"""
    suitability(layers, index, t::Number)
Returns a scalar representing cell suitability from one or multiple layers.

For multiple layers the product of each individual scalar is returned.

Layers of type other than AbstractLayer return 1.0.

### Arguments
- `layers` : a single layer or tuple of layers of any type
- `index` : a tuple containing 2 grid coordinates
- `t::Number` : current timestep for interploating layere sequences.
"""
function get_layers end
get_layers(model, data, index) = get_layers(model, data, model.layers, index, data.t)
get_layers(model, data, layers::Tuple, index, t::Number) =
    get_layers(model, data, layers[1], index, t) *
    get_layer(model, data, Base.tail(layers), index, t::Number)
get_layers(model, data, layers::Tuple{}, index, t::Number) = 1
get_layers(model, data, layer::Matrix, index, t::Number) =
    return layer[index...]
get_layers(model, data, sequence::AbstractSequence, index, t::Number) =
    sequence_interpolate(model, data, sequence, index, t)
get_layers(model, data, layers, index, t::Number) = 1

"""
    sequence_interpolate(layer, index, t)
Interpolates between layers in a sequence.
"""
@inline sequence_interpolate(model, data, seq, index, t::Number) = begin
    # Linear interpolation between layers in the sequence.
    tf = calc_timestep(data, seq, t)
    tr = unsafe_trunc(Int64, tf)
    frac = tf - tr
    pos1 = cyclic(tr, length(seq))
    pos2 = cyclic(tr + oneunit(tr), length(seq))
    seq[pos1][index...] * (oneunit(frac) - frac) + seq[pos2][index...] * frac
end

# Time position is centered in the current frame, not at the start.
calc_timestep(data, seq, t) = t * (data.timestep / seq.timestep) + oneunit(t) / 2

"Cycles a time position through a particular timestep length"
cyclic(t::Int, len::Int) =
    if t > len
        rem(t + len - 1, len) + 1
    elseif t <= 0
        t + (t ÷ len -1) * -len
    else
        t
    end

abstract type Layer end
abstract type Sequence <: Layer end

struct SuitabilityLayer <: Layer end
struct SuitabilitySequence <: Sequence end

struct HumanLayer <: Layer end
struct HumanSequence <: Sequence end

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
suitability(layers::Tuple, index, t::Number) = mapreduce(l -> suitability(l, index, t), *, layers)
suitability(layer, index, t::Number) = 1.0
suitability((X, layer)::Tuple{SuitabilityLayer, AbstractMatrix}, index, t::Number) = layer[index...]
suitability((_, sequence, timespan)::Tuple{SuitabilitySequence, Tuple, Number}, index, t::Number) = 
    sequence_interpolate(sequence, timespan, index, t)


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
    mapreduce(l -> human_impact(l, index, t), *, layers)
human_impact(::Tuple{Layer,Vararg}, index, t::Number) = 1.0
human_impact((_, layer)::Tuple{HumanLayer,Tuple}, index, t) = layer[index...]

"""
    sequence_interpolate(layer, index, t)
Interpolates between layers in a sequence.
"""
sequence_interpolate(sequence, timespan, index, t::Number) = begin
    # Time position is centered in the current frame, not at the start.
    tf = (t + timespan * 0.5) / timespan
    # Linear interpolation between layers in the sequence.
    tr = floor(Int64, tf)
    frac = convert(eltype(sequence[1]), tf - tr)
    pos1 = cyclic(tr, length(sequence))
    pos2 = cyclic(tr + 1, length(sequence))
    @inbounds sequence[pos1][index...] * (oneunit(eltype(sequence[1])) - frac) + sequence[pos2][index...] * frac
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

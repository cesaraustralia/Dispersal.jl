# TODO a better way to do this.
# Probably using traits. Which layers are included should be flixible
# - you should be able to pass in any struct with the right field names and it
# will do all the right things. Single inheritance can't capture that.
"""
Layers wrap an array, normally built from some kind of raster.

The wrapper defines the purpose of the layer and allows specialised
method dispatch to utilise them.
"""
abstract type AbstractLayer end

abstract type AbstractHumanLayer <: AbstractLayer end
"""
A wrapper for arrays that provide human dispersal scalars for the grid.
"""
struct HumanLayer{T} <: AbstractHumanLayer
    "Any 2-dimensional AbstractArray matching the coordinates of the init array"
    data::T
end

"""
Abstract type for layers that provide suitability scalars.
"""
abstract type AbstractSuitabilityLayer <: AbstractLayer end

"""
A wrapper for arrays that provide suitability scalars for the grid.
"""
struct SuitabilityLayer{T} <: AbstractSuitabilityLayer
    "Any 2-dimensional AbstractArray matching the coordinates of the init array"
    data::T
end

"""
An abstract type for layers that provide sequences of suitability scalars.
"""
abstract type AbstractSuitabilitySequence <: AbstractSuitabilityLayer end

"""
A wrapper for a sequence of arrays that provide suitability scalar
for grid points that cycle over a time span.
"""
struct SuitabilitySequence{T,D} <: AbstractSuitabilitySequence
    "The timespan of each layer in the sequence."
    timespan::T
    """
    Either an Array of 2-dimensional arrays matching the coordinates of the init array,
    or a similar 3 dimensional array where the 3rd dimension is the time-step.
    """
    data::D
end


"""
    suitability(layers, row::Int, col::Int, t::Number)
Returns a scalar representing cell suitability from one or multiple layers.

For multiple layers the product of each individual scalar is returned.

Layers of type other than AbstractSuitabilityLayer return 1.0.

### Arguments
- `layers` : a single layer or tuple of layers of any type
- `row::Int`
- `col::Int`
- `t::Number` : current timestep for interploating layere sequences.
"""
function suitability end
suitability(layers, row::Int, col::Int, t::Number) = 1.0
suitability(layers::Tuple, row::Int, col::Int, t::Number) =
    mapreduce(l -> suitability(l, row, col, t), *, 1.0, layers)
suitability(layer::AbstractSuitabilityLayer, row::Int, col::Int, t::Number) =
    get_cell(layer, row, col, round(Int,t))
suitability(layer::AbstractSuitabilitySequence, row::Int, col::Int, t::Number) =
    sequence_interpolate(layer, row, col, t)


"""
    human_impact(layers, row::Int, col::Int, t)
Returns a scalar indicating human impact from one or multiple layers.

For multiple layers the product of each individual scalar is returned.

Layers of type other than AbstractHumanLayer return 1.0.

### Arguments
- `layers` : a single layer or tuple of layers of any type
- `row::Int`
- `col::Int`
- `t::Number` : current timestep for interploating layere sequences.
"""
function human_impact end
human_impact(layers, row::Int, col::Int, t) = 1.0
human_impact(layers::Tuple, row::Int, col::Int, t) =
    mapreduce(f -> human_impact(getfield(layers, f), row, col, t), *, 1.0, fieldnames(L))
human_impact(layer::AbstractHumanLayer, row::Int, col::Int, t) =
    get_cell(layer, row, col, round(Int,t))

"""
    get_cell(layer, row::Int, col::Int, pos::Number)
Return a particular cell from a layer, given row, column and timestep)
"""
get_cell(layer, row::Int, col::Int, t::Number) = get_cell(layer.data, row::Int, col::Int, t)
get_cell(data::AbstractArray{A,1}, row::Int, col::Int, t::Number) where A <: AbstractArray = data[t][row, col]
get_cell(data::AbstractArray{T,2}, row::Int, col::Int, t::Number) where T = data[row, col]
get_cell(data::AbstractArray{T,3}, row::Int, col::Int, t::Number) where T = data[row, col, t]

timespan(layer) = layer.timespan

"""
    sequence_interpolate(layer, row, col, t)
Interpolates between layers in a sequence.
"""
sequence_interpolate(layer, row::Int, col::Int, t::Number) = begin
    # Time position is centered in the current frame, not at the start.
    tf = (t + timespan(layer) * 0.5) / timespan(layer)
    # Linear interpolation between layers in the sequence.
    tr = floor(Int, tf)
    frac = tf - tr
    pos1 = cyclic(tr, num_frames(layer))
    pos2 = cyclic(tr + 1, num_frames(layer))
    get_cell(layer, row, col, pos1) * (1.0 - frac) + get_cell(layer, row, col, pos2) * frac
end

num_frames(layer::AbstractLayer) = num_frames(layer.data)
num_frames(data::AbstractArray{A,1}) where A <: AbstractArray = length(data)
num_frames(data::AbstractArray{T,3}) where T = size(data, 3)

"Cycles a time position through a particular timespan length"
cyclic(t::Number, len::Number) =
    if t > len
        t - (t ÷ len) * len
    elseif t <= 0
        t + (t ÷ len -1) * -len
    else
        t
    end


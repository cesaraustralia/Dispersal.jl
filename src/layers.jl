# TODO a better way to do this.
# Probably using traits. Which layers are included should be flixible
# - you should be able to pass in any struct with the right field names and it
# will do all the right things. Single inheritance can't capture that.
"""
$(TYPEDEF)
Layers wrap an array, normally built from some kind of raster.

The wrapper defines the purpose of the layer and allows specialised 
method dispatch to utilise them.
"""
abstract type AbstractLayer end

"""
$(TYPEDEF)
$(FIELDS)
"""
struct HumanLayer{T} <: AbstractLayer
    data::T
end

"""
$(TYPEDEF)
"""
abstract type AbstractSuitabilityLayer <: AbstractLayer end
"""
$(TYPEDEF)
$(FIELDS)
"""
struct SuitabilityLayer{T} <: AbstractSuitabilityLayer
    data::T
end

"""
$(TYPEDEF)
"""
abstract type AbstractSuitabilitySequence <: AbstractSuitabilityLayer end

"""
$(TYPEDEF)
$(FIELDS)
"""
struct SuitabilitySequence{T,D} <: AbstractSuitabilitySequence
    timespan::T
    data::D
end

"""
$(TYPEDEF)
"""
abstract type AbstractLayers end

"""
$(TYPEDEF)
$(FIELDS)
"""
struct Layers{S,F} <: AbstractLayers
    suitability::S
    human::F
end


"""
$(SIGNATURES)
Returns a suitability scalar from a single layer, or the product of multiple layers
"""
suitability(layers::L, row::Int, col::Int, t) where L <: AbstractLayers =
    mapreduce(f -> suitability(getfield(layers, f), row, col, t), *, 1.0, fieldnames(L))
suitability(layer::AbstractSuitabilityLayer, row::Int, col::Int, t) = 
    get_cell(layer, row, col, round(Int,t))
suitability(layer::AbstractSuitabilitySequence, row::Int, col::Int, t) = 
    seq_interpolate(layer, row, col, t)
suitability(layers, row::Int, col::Int, t) = 1.0

frames(layer::AbstractLayer) = frames(layer.data)
frames(data::AbstractArray{A,1}) where A <: AbstractArray = length(data)
frames(data::AbstractArray{T,3}) where T = size(data, 3)

"""
$(SIGNATURES)
Return a particular cell from a layer, given row, column and timestep)
"""
get_cell(layer, row::Int, col::Int, pos) = get_cell(layer.data, row::Int, col::Int, pos)
get_cell(data::AbstractArray{A,1}, row::Int, col::Int, pos) where A <: AbstractArray = data[pos][row, col]
get_cell(data::AbstractArray{T,2}, row::Int, col::Int, pos) where T = data[row, col]
get_cell(data::AbstractArray{T,3}, row::Int, col::Int, pos) where T = data[row, col, pos]

timespan(layer) = layer.timespan

"""
$(SIGNATURES)
Interpolates between layers in a sequence
"""
seq_interpolate(layer, row, col, t) = begin
    tf = (t + timespan(layer) * 0.5) / timespan(layer) # position centered in frame
    # linear interpolation between layers in the sequence
    tr = floor(Int, tf)
    frac = tf - tr
    pos1 = cyclic(tr, frames(layer))
    pos2 = cyclic(tr + 1, frames(layer))
    get_cell(layer, row, col, pos1) * (1.0 - frac) + get_cell(layer, row, col, pos2) * frac
end

cyclic(t::Int, len::Int) = 
    if t > len
        1
    elseif t <= 0
        len
    else
        t
    end



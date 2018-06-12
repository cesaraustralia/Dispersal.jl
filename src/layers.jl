# TODO a better way to do this.
# Probably using traits. Which layers are included should be flixible
# - you should be able to pass in any struct with the right field names and it
# will do all the right things. Single inheritance can't capture that.
abstract type AbstractLayer end
struct HumanLayer{T} <: AbstractLayer
    data::T
end
abstract type AbstractSuitabilityLayer <: AbstractLayer end
struct SuitabilityLayer{T} <: AbstractSuitabilityLayer
    data::T
end

abstract type AbstractSuitabilitySequence <: AbstractSuitabilityLayer end
struct SuitabilitySequence{T,D} <: AbstractSuitabilitySequence
    timespan::T
    data::D
end

abstract type AbstractLayers end
struct Layers{S,F} <: AbstractLayers
    suitability::S
    human::F
end


"""
Returns a suitability scalar from a single layer, or the product of multiple layers
"""
suitability(layers::L, row::Int, col::Int, t) where L <: AbstractLayers =
    mapreduce(f -> suitability(getfield(layers, f), row, col, t), *, 1.0, fieldnames(L))

suitability(layer::AbstractSuitabilityLayer, row::Int, col::Int, t) = layer_coord(layer.data, row, col, round(Int,t))
suitability(layer::AbstractSuitabilitySequence, row::Int, col::Int, t) = begin
    tf = (t + layer.timespan * 0.5) / layer.timespan # position centered in frame
    # linear interpolation between layers in the sequence
    tr = floor(Int, tf)
    frac = tf - tr
    pos1 = cyclic(tr, frames(layer.data))
    pos2 = cyclic(tr + 1, frames(layer.data))
    layer_coord(layer.data, row, col, pos1) * (1.0 - frac) + layer_coord(layer.data, row, col, pos2) * frac
end
suitability(layers, row::Int, col::Int, t) = 1.0

frames(data::AbstractArray{A,1}) where A <: AbstractArray = length(data)
frames(data::AbstractArray{T,3}) where T = size(data, 3)

layer_coord(data::AbstractArray{A,1}, row::Int, col::Int, pos) where A <: AbstractArray = data[pos][row, col]
layer_coord(data::AbstractArray{T,2}, row::Int, col::Int, pos) where T = data[row, col]
layer_coord(data::AbstractArray{T,3}, row::Int, col::Int, pos) where T = data[row, col, pos]

cyclic(t::Int, len::Int) = 
    if t > len
        1
    elseif t <= 0
        len
    else
        t
    end



@mix @columns struct Layers{L,TI}
    layer::L           | nothing | false | _          | "Data layer"
    timeindex::TI      | 1       | false | _          | "Precalculated interpolation indices"
end

"""
    layer(rule)

Returns the value of a single layer or interplated value from a sequence of layers.

If multiple layers are available the product will be returned.
"""
function layer end
layer(rule::Rule) = rule.layer
Base.@propagate_inbounds layer(rule::Rule, data, index) =
    layer(layer(rule), data, index, rule.timeindex)
Base.@propagate_inbounds layer(l::Matrix, data, index, timeindex) = l[index...]
Base.@propagate_inbounds layer(l::AbstractArray{T,3}, data, index, timeindex) where T =
    l[index..., timeindex]

# Layer precalcs
precalclayer(::AbstractMatrix, rule::Rule, data) = rule
precalclayer(::AbstractArray{<:Any,3}, rule::Rule, data) =
    @set rule.timeindex = precalc_timeindex(layer(rule), rule, data)

"""
    precalc_timeindex(layer, index, t)

Interpolates between layers in a sequence. This should probably be
replaced with an external interpolation package.
"""
# Truncated version
precalc_timeindex(layer, rule, data, t=currenttime(data)) = begin
    # Convert Month etc timesteps to a realised DateTime period
    layerstep = length(starttime(layer):timestep(layer):t)
    cyclic_index(layerstep, size(layer, Time))
end
# Interpolated version
# Base.@propagate_inbounds precalc_timeindex(layers, rule, data) = begin
#     # Convert Month etc timesteps to a realised DateTime period
#     f = currentframe(data)
#     sim_timestep = t - (t + timestep(data))
#     layer_timestep = t - (t + timestep(layers))
#     # Time position is centered in the current frame, not at the start.
#     # this allows interpolating between two positions.
#     # Needs multiple oneunit() to avoid promotion to Float64
#     t_float = t * (sim_timestep / layer_timestep) + 0.5
#     t_int = unsafe_trunc(Int64, t_float)
#     frac = t_float - t_int
#     len = size(layers, Time)
#     t1 = cyclic(t_int, len)
#     t2 = cyclic(t_int + oneunit(t_int), len)
#     WeightedArbIndex((t1, t2), (frac, 1-frac))
# end

"Cycles a time position through a particular timestep length"
cyclic_index(i::Integer, len::Integer) =
    if i > len
        rem(i + len - 1, len) + 1
    elseif i <= 0
        i + (i ÷ len -1) * -len
    else
        i
    end


# Get time step/start/stop from AbstractDimensionalArray
# Integer fallbacks are for other array types is the indices of dim 3
DynamicGrids.timestep(A::AbstractDimensionalArray) = step(dims(A, Time))
DynamicGrids.timestep(A::AbstractArray) = 1

DynamicGrids.starttime(A::AbstractDimensionalArray) = first(dims(A, Time))
DynamicGrids.starttime(A::AbstractArray) = firstindex(A, 3)

DynamicGrids.stoptime(A::AbstractDimensionalArray) = last(dims(A, Time))
DynamicGrids.stoptime(A::AbstractArray) = lastindex(A, 3)


"""
Simple rule that copies a layer to a grid over time. This can be used for
comparing simulation dynamics to layer dynamics.
"""
@Layers struct LayerCopy{L} <: Rule end

DynamicGrids.applyrule(rule::LayerCopy, data, state, index, args...) =
    layer(rule, data, index)

DynamicGrids.precalcrules(rule::LayerCopy, data) =
    precalclayer(layer(rule), rule, data)

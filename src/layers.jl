
@mix @columns struct Layers{L,TI}
    layerkey::L   | nothing | false | _ | "key for aux layer"
    timeindex::TI | 1       | false | _ | "Precalculated interpolation indices. Not set by users"
end

"""
layer(rule::Rule, data::SimData, [I...])

Returns the value of a single layer or interplated value from a sequence of layers.

If multiple layers are available the product will be returned. Corresponding
layers must be include in the `aux` `NamedTuple` in the out put or passed to `sim!`
"""
function layer end
layer(rule::Rule, data) = aux(data)[unwrap(layerkey(rule))]
layer(rule::Rule, data, I) =
    layer(layer(rule, data), I, timeindex(rule))
layer(rule::Rule, data, I, L) =
    layer(layer(rule, data), I, timeindex(rule), L)

layer(l::Matrix, I, timeindex) = l[I...]
layer(l::AbstractArray{T,3}, I, timeindex) where T =
    l[I..., timeindex]
layer(l::AbstractArray{T,4}, I, timeindex, layerindex) where T =
    l[I..., timeindex, layerindex]

timeindex(rule::Rule) = rule.timeindex
layerkey(rule::Rule) = rule.layerkey

# Layer precalcs
precalclayer(::AbstractMatrix, rule::Rule, data) = rule
precalclayer(::AbstractArray{<:Any,3}, rule::Rule, data) =
    @set rule.timeindex = precalc_timeindex(layer(rule, data), rule, data)
precalclayer(::AbstractArray{<:Any,4}, rule::Rule, data) =
    @set rule.timeindex = precalc_timeindex(layer(rule, data), rule, data)
"""
    precalc_timeindex(layer, index, t)

Interpolates between layers in a sequence. This should probably be
replaced with an external interpolation package.
"""
# Truncated version
precalc_timeindex(layer, rule, data, t=currenttime(data)) = begin
    # Convert Month etc timesteps to a realised DateTime period
    layerstep = length(starttime(layer):timestep(layer):t)
    cyclic_index(layerstep, size(layer, 3))
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
#     t_float = t * (sim_timestep / layer_timestep) + 0.5 #     t_int = unsafe_trunc(Int64, t_float)
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
DynamicGrids.timestep(A::AbstractDimensionalArray) = step(dims(A, TimeDim))
DynamicGrids.timestep(A::AbstractArray) = 1

starttime(A::AbstractDimensionalArray) = first(dims(A, TimeDim))
starttime(A::AbstractArray) = firstindex(A, 3)

stoptime(A::AbstractDimensionalArray) = last(dims(A, TimeDim))
stoptime(A::AbstractArray) = lastindex(A, 3)

"""
    LayerCopy(layerkey, timeindex)

A simple rule that copies a layer to a grid over time. 
This can be used for comparing simulation dynamics to layer dynamics.

$(FIELDDOCTABLE)
"""
@Layers struct LayerCopy{R,W} <: Rule{R,W} end
LayerCopy(args...) = LayerCopy{Tuple{},:_default_,map(typeof,args)...}(args...)

DynamicGrids.applyrule(data, rule::LayerCopy, state, index, layerindex, args...) =
    layer(rule, data, index, layerindex)

DynamicGrids.applyrule(data, rule::LayerCopy, state, index, args...) =
    layer(rule, data, index)

DynamicGrids.precalcrules(rule::LayerCopy, data) =
    precalclayer(layer(rule, data), rule, data)

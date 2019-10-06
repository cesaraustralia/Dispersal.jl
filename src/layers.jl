
timestep(a::AbstractDimensionalData) = step(val(dims(a, Time)))
starttime(a::AbstractDimensionalData) = first(val(dims(a, Time)))
stoptime(a::AbstractDimensionalData) = last(val(dims(a, Time)))

"""
    layer(rule)

Returns the value of a single layer or interplated value from a sequence of layers.

If multiple layers are available the product will be returned.
"""
function layer end
layer(rule::AbstractRule) = rule.layer
Base.@propagate_inbounds layer(rule::AbstractRule, data, index) =
    layer(layer(rule), data, index, timeinterp(rule))
Base.@propagate_inbounds layer(layer::Matrix, data, index, interp) = layer[index...]
Base.@propagate_inbounds layer(layers::AbstractArray{T,3}, data, index, interp) where T =
    layers[index..., interp]

"""
    precalc_time_interpolation(layer, index, t)

Interpolates between layers in a sequence. This should probably be 
replaced with an external interpolation package.
"""
precalc_time_interpolation(layer::AbstractMatrix, rule, data) = currentframe(data) 
# Base.@propagate_inbounds precalc_time_interpolation(layers, rule, data) = begin
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
Base.@propagate_inbounds precalc_time_interpolation(layers, rule, data, t = currenttime(data)) = begin
    # Convert Month etc timesteps to a realised DateTime period
    tstep = t + timestep(layers) - t
    nsteps = length(starttime(layer):tstep:t)
    len = size(layers, Time)
    cyclic(t_int, len)
end

"Cycles a time position through a particular timestep length"
cyclic(i::Integer, len::Integer) =
    if i > len
        rem(i + len - 1, len) + 1
    elseif i <= 0
        i + (i ÷ len -1) * -len
    else
        i
    end

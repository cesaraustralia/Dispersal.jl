timestep(a) = step(val(dims(a, Time)))

"""
    layer(rule)

Returns the value of a single layer or interplated value from a sequence of layers.

If multiple layers are available the product will be returned.
"""
function layer end
layer(rule::AbstractRule) = rule.layer
Base.@propagate_inbounds layer(rule::AbstractRule, data, index) =
    layer(layer(rule), data, index, timeinterp(rule))
Base.@propagate_inbounds layer(layer::Matrix, data, index, t) = layer[index...]
Base.@propagate_inbounds layer(layers::AbstractArray{T,3}, data, index, t) where T =
    interpolate_layers(layers, index, t)

"""
    interpolate_layers(layers, index, t)

Interpolates between layers in a sequence.
"""
Base.@propagate_inbounds interpolate_layers(layers, index, (t1, t2, frac)) = 
    layers[index..., t1] * (oneunit(frac) - frac) + layers[index..., t2] * frac

"""
    precalc_time_interpolation(layer, index, t)

Interpolates between layers in a sequence. This should probably be 
replaced with an external interpolation package.
"""
precalc_time_interpolation(layer::Matrix, rule, data) = timeinterp(rule)
Base.@propagate_inbounds precalc_time_interpolation(layers, rule, data) = 
    precalc_time_interpolation(layers, rule, data, timestep(data), currenttime(data))
Base.@propagate_inbounds precalc_time_interpolation(layers, rule, data, sim_timestep, t) = begin
    # Time position is centered in the current frame, not at the start.
    # this allows interpolating between two positions.
    layer_timestep = timestep(layers)
    # Needs multiple oneunit() to avoid promotion to Float64
    tf = t * (sim_timestep / layer_timestep) + oneunit(t) / 2one(sim_timestep)
    tr = unsafe_trunc(Int64, tf)
    frac = tf - tr
    len = size(layers, Time)
    t1 = cyclic(tr, len)
    t2 = cyclic(tr + oneunit(tr), len)
    t1, t2, frac
end

"Cycles a time position through a particular timestep length"
cyclic(t::Integer, len::Integer) =
    if t > len
        rem(t + len - 1, len) + 1
    elseif t <= 0
        t + (t ÷ len -1) * -len
    else
        t
    end

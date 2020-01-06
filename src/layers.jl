
DynamicGrids.timestep(A::AbstractDimensionalArray) = step(dims(A, Time))
DynamicGrids.starttime(A::AbstractDimensionalArray) = first(val(dims(A, Time)))
DynamicGrids.stoptime(A::AbstractDimensionalArray) = last(val(dims(A, Time)))
# Interger gallback for other array types is the indices of dim 3
DynamicGrids.timestep(A::AbstractArray) = 1
DynamicGrids.starttime(A::AbstractArray) = firstindex(A, 3)
DynamicGrids.stoptime(A::AbstractArray) = lastindex(A, 3)

"""
    layer(rule)

Returns the value of a single layer or interplated value from a sequence of layers.

If multiple layers are available the product will be returned.
"""
function layer end
layer(rule::Rule) = rule.layer
Base.@propagate_inbounds layer(rule::Rule, data, index) =
    layer(layer(rule), data, index, rule.timeinterp)
Base.@propagate_inbounds layer(l::Matrix, data, index, interp) = l[index...]
Base.@propagate_inbounds layer(l::AbstractArray{T,3}, data, index, interp) where T =
    l[index..., interp]

# Layer precalcs
precalclayer(::AbstractMatrix, rule::Rule, data) = rule
precalclayer(::AbstractArray{<:Any,3}, rule::Rule, data) =
    @set rule.timeinterp = precalc_time_interpolation(layer(rule), rule, data)

"""
    precalc_time_interpolation(layer, index, t)

Interpolates between layers in a sequence. This should probably be
replaced with an external interpolation package.
"""
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
Base.@propagate_inbounds precalc_time_interpolation(layer, rule, data, t=currenttime(data)) = begin
    # Convert Month etc timesteps to a realised DateTime period
    tstep = t + timestep(layer) - t
    nsteps = length(starttime(layer):tstep:t)
    len = size(layer, Time)
    cyclic(nsteps, len)
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


precalc_timestep(rule, data) = precalc_timestep(rule.timestep, rule, data)
precalc_timestep(ruletimestep::DatePeriod, rule, data) =
    @set rule.nsteps = currenttimestep(data) / Millisecond(ruletimestep)
precalc_timestep(ruletimestep::Nothing, rule, data) = @set rule.nsteps = 1
precalc_timestep(ruletimestep, rule, data) =
    @set rule.nsteps = timestep(data) / ruletimestep

"""
    precalc_auxtimeindex(A, rule, data)

Sets the auxtimeindex field based on the auxilary data time axis, 
and the current simulation timestep.

Returns the updated rule, or an unchanged rule if the aux data
ia a Matrix, without a time dimension.
"""
function precalc_auxtimeindex(
    A::AbstractArray{<:Any,3}, rule::Rule, data::SimData, t=currenttime(data)
)
    # Convert Month etc timesteps to a realised DateTime period
    step = length(starttime(A):timestep(A):t)
    @set rule.auxtimeindex = cyclic_index(step, size(A, 3))
end
precalc_auxtimeindex(::AbstractMatrix, rule::Rule, data::SimData, t...) = rule


"""
    auxval(rule::Rule, data::SimData, key, I...)

Returns the value of a single layer or interplated value from a sequence of layers.

Corresponding layers must be include as the `aux` keyword to the `Output` or `sim!`.
"""
function auxval end
auxval(data::SimData, key::Union{Symbol,Val}, I...) = auxval(aux(data, key), I...) 
# If there is no time dimension we return the same data for every timestep
auxval(A::Matrix, y, x, t=nothing) = A[y, x] 
# If there is a time dimension we index into it with `t`
auxval(A::AbstractArray{T,3}, y, x, t) where T = A[y, x, t]

# Get time step/start/stop from AbstractrimArray
# For other array types the fallback is the indices of dim 3
DynamicGrids.timestep(A::AbstractDimArray) = step(dims(A, TimeDim))
DynamicGrids.timestep(A::AbstractArray) = 1

starttime(A::AbstractDimArray) = first(dims(A, TimeDim))
starttime(A::AbstractArray) = firstindex(A, 3)

stoptime(A::AbstractDimArray) = last(dims(A, TimeDim))
stoptime(A::AbstractArray) = lastindex(A, 3)

# Interpolated version
# Base.@propagate_inbounds precalc_auxtimeindex(A, rule, data) = begin
#     # Convert Month etc timesteps to a realised DateTime period
#     f = currentframe(data)
#     sim_timestep = t - (t + timestep(data))
#     layer_timestep = t - (t + timestep(A))
#     # Time position is centered in the current frame, not at the start.
#     # this allows interpolating between two positions.
#     # Needs multiple oneunit() to avoid promotion to Float64
#     t_float = t * (sim_timestep / layer_timestep) + 0.5 #     t_int = unsafe_trunc(Int64, t_float)
#     frac = t_float - t_int
#     len = size(A, Time)
#     t1 = cyclic(t_int, len)
#     t2 = cyclic(t_int + oneunit(t_int), len)
#     WeightedArbIndex((t1, t2), (frac, 1-frac))
# end

"Cycles a time position through a particular timestep length"
function cyclic_index(i::Integer, len::Integer)
    return if i > len
        rem(i + len - 1, len) + 1
    elseif i <= 0
        i + (i ÷ len -1) * -len
    else
        i
    end
end

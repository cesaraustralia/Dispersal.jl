"Extend to modify [`JumpDispersal`](@ref)"
abstract type AbstractJumpDispersal <: AbstractPartialModel end

"Jump dispersal within a [`DispersalNeighborhood`](@ref)] or other neighborhoods."
@Probabilistic struct JumpDispersal{SR} <: AbstractJumpDispersal 
    "A number or Unitful.jl distance with the same units as cellsize"
    spotrange::SR = 30.0 | true | (0.0, 100.0)
end

"""
    rule(model::AbstractJumpDispersal, data, state, index, args...)
Long range rule for [`AbstractJumpDispersal`](@ref). A random cell
within the spotrange is invaded if it is suitable.
"""
@inline rule!(model::AbstractJumpDispersal, data, state, index, args...) = begin
    # Ignore empty cells
    state > zero(state) || return state

    # Random dispersal events
    spec_rand(data.source, Float64, args...) < model.prob_threshold || return state

    # Randomly select spotting distance
    rnge = spec_rand.((data.source,), (Float64, Float64), tuple.(args)...) .* (model.spotrange / data.cellsize)
    spot = tuple(unsafe_trunc.(Int64, rnge .+ index)...)
    spot, is_inbounds = inbounds(spot, data.dims, Skip())

    # Update spotted cell if it's on the grid
    if is_inbounds
        @inbounds data.dest[spot...] = state
    end

    state
end

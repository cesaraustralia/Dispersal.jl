"Extends AbstractPartialRule for jump dispersal rules"
abstract type AbstractJumpDispersal <: AbstractPartialRule end

"""
Jump dispersal simulates random long distance dispersal events. A random cell within 
the spotrange is invaded if it is suitable.  A [`Mask`](@ref) rule may be useful 
aftwer this rule, as dispersal events may be to anywhere on the grid within the given range.
$(FIELDDOCTABLE)
"""
@Probabilistic struct JumpDispersal{SR} <: AbstractJumpDispersal 
    # Field       | Def  | Flatten | Limits       | Description
    spotrange::SR | 30.0 | true    | (0.0, 100.0) | "A number or Unitful.jl distance with the same units as cellsize"
end

@inline applyrule!(rule::AbstractJumpDispersal, data, state, index, args...) = begin
    # Ignore empty cells
    state > zero(state) || return state

    # Random dispersal events
    rand() < rule.prob_threshold || return state

    # Randomly select spotting distance
    rnge = rand(2) .* (rule.spotrange / data.cellsize)
    spot = tuple(unsafe_trunc.(Int64, rnge .+ index)...)
    spot, is_inbounds = inbounds(spot, data.dims, Skip())

    # Update spotted cell if it's on the grid
    if is_inbounds
        @inbounds data.dest[spot...] = state
    end

    state
end

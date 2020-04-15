"""
Abstract supertype that extends [`PartialRule`](@ref) for jump dispersal rules
"""
abstract type AbstractJumpDispersal{R,W} <: PartialRule{R,W} end

"""
    JumpDispersal(spotrange)
    JumpDispersal{R,W}(spotrange)
    JumpDispersal(; spotrange=30.0)

Jump dispersal simulates random long distance dispersal events. A random cell within 
the `spotrange` is invaded. 

Pass grid name `Symbol`s to `R` and `W` type parameters to use specific grids.

$(FIELDDOCTABLE)
"""
@Probabilistic struct JumpDispersal{R,W,SR} <: AbstractJumpDispersal{R,W}
    # Field       | Def  | Flatten | Limits       | Description
    spotrange::SR | 30.0 | true    | (0.0, 100.0) | "A number or Unitful.jl distance with the same units as cellsize"
end

# TODO update this and test
@inline applyrule!(rule::AbstractJumpDispersal{R,W}, data, state, index
                  ) where {R,W} = begin
    # Ignore empty cells
    state > zero(state) || return state

    # Random dispersal events
    rand() < rule.prob_threshold || return state

    # Randomly select spotting distance
    rnge = rand(2) .* rule.spotrange
    spot = tuple(unsafe_trunc.(Int64, rnge .+ index)...)
    spot, is_inbounds = inbounds(spot, gridsize(data), RemoveOverflow())

    # Update spotted cell if it's on the grid
    if is_inbounds
        @inbounds dest(data)[W][spot...] = state
    end

    state
end

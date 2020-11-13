"""
Jump dispersal rules.
"""
abstract type AbstractJumpDispersal{R,W} <: ManualRule{R,W} end

"""
    JumpDispersal(spotrange)
    JumpDispersal{R,W}(spotrange)
    JumpDispersal(; spotrange=30.0)

Jump dispersal simulates random long distance dispersal events. A random cell within 
the `spotrange` is invaded. 

Pass grid name `Symbol`s to `R` and `W` type parameters to use specific grids.
"""
struct JumpDispersal{R,W,PT,SR} <: AbstractJumpDispersal{R,W}
    "A real number between one and zero"
    prob_threshold::PT
    "A number or Unitful.jl distance with the same units as cellsize"
    spotrange::SR
end

JumpDispersal{R,W}(; 
    prob_threshold=Param(0.1, bounds=(0.0, 1.0)),
    spotrange=Param(30.0, bounds=(0.0, 100.0)),
    ) where {R,W} = JumpDispersal{R,W}(prob_threshold, spotrange)

# DynamicGrids.jl interface

# TODO update this and test
@inline function applyrule!(data, rule::AbstractJumpDispersal{R,W},
                            state, cellindex) where {R,W}
    # Ignore empty cells
    state > zero(state) || return state

    # Random dispersal events
    rand() < rule.prob_threshold || return state

    # Randomly select spotting distance
    intspot = round(Int, rule.spotrange)
    rnge = -intspot:intspot
    jump = (rand(rnge), rand(rnge))
    jumpdest, is_inbounds = inbounds(jump .+ cellindex, gridsize(data), RemoveOverflow())

    # Update spotted cell if it's on the grid
    is_inbounds && @inbounds add!(data[W], state, jumpdest...)

    return state
end

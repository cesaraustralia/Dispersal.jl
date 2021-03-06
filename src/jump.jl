"""
    JumpDispersal <: SetCellRule

    JumpDispersal(; prob_threshold, spotrange)
    JumpDispersal{R}(; prob_threshold, spotrange)
    JumpDispersal{R,W}(; prob_threshold, spotrange)

Jump dispersal simulates a random long distance dispersal event, which occurs at a given
probablity. A random cell within the `spotrange` is invaded by duplicating the entire 
population of the source cell at the destination cell. 

# Keyword Arguments

- `prob_threshold`: probability of dispersal, between one and zero. 
    May be a `Number`, an `Aux` array or another `Grid`:
- `spotrange`: number of cells in range of jumps, in any direction.

Pass grid `Symbol`s to `R` or both `R` and `W` type parameters to use to specific grids.
"""
struct JumpDispersal{R,W,PT,SR} <: SetCellRule{R,W}
    prob_threshold::PT
    spotrange::SR
end
function JumpDispersal{R,W}(; 
    prob_threshold=Param(0.1, bounds=(0.0, 1.0)),
    spotrange=Param(30.0, bounds=(0.0, 100.0)),
) where {R,W}
    JumpDispersal{R,W}(prob_threshold, spotrange)
end

@inline function applyrule!(data, rule::JumpDispersal{R,W}, N, I) where {R,W}
    # Ignore empty cells
    N > zero(N) || return nothing
    # Random dispersal events
    p = get(data, rule.prob_threshold, I...) 
    rand() < p || return nothing

    # Randomly select spotting distance
    intspot = round(Int, rule.spotrange)
    rnge = -intspot:intspot
    dest, is_inbounds = inbounds(data, (rand(rnge), rand(rnge)) .+ I)

    # Update spotted cell if it's on the grid
    is_inbounds && @inbounds add!(data[W], N, dest...)

    return nothing
end

"""
    OutwardsDispersal(neighborhood)
    OutwardsDispersal(; neighborhood=DispersalKernel{3}())
    OutwardsDispersal{R,W}(neighborhood)
    OutwardsDispersal{R,W}(; neighborhood=DispersalKernel{3}())

Dispersal reduces the current cell object, increasing the populations of the
cells in the surrounding neighborhood.

Pass grid name `Symbol`s to `R` and `W` type parameters to use specific grids.

The result should be identical to [`InwardsDispersal`](@ref) but may be more
efficient than when a small number of cells are occupied. It is less efficient
when a large proportion of the grid is occupied.
"""
struct OutwardsDispersal{R,W,N} <: SetNeighborhoodRule{R,W}
    neighborhood::N
end
function OutwardsDispersal{R,W}(; neighborhood=DispersalKernel{3}()) where {R,W} 
    OutwardsDispersal{R,W}(neighborhood)
end

@inline function applyrule!(data, rule::OutwardsDispersal{R,W}, state, I) where {R,W}
    state == zero(state) && return
    hood = neighborhood(rule)
    sum = zero(state)
    for offset in offsets(rule)
        @inbounds propagules = state * kernel(hood)[(offset .+ radius(hood) .+ 1)...]
        @inbounds add!(data[W], propagules, I .+ offset...)
        sum += propagules
    end
    # Subtract from current cell, unless state is Bool
    state isa Bool || @inbounds sub!(data[W], sum, I...)
    return nothing
end

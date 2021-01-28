"""
    OutwardsPopulationDispersal(neighborhood)
    OutwardsPopulationDispersal(; neighborhood)
    OutwardsPopulationDispersal{R}(; neighborhood)
    OutwardsPopulationDispersal{R,W}(; neighborhood)

Dispersal reduces the current cell population, increasing the populations of the
cells in the surrounding neighborhood deterministically.

This will make sense ecologically where cell populations are large, 
otherwise a randomised kernel may be more suitable.

The result should be identical to [`InwardsDispersal`](@ref) but may be more
efficient than when a small number of cells are occupied. It is less efficient
when a large proportion of the grid is occupied.

# Arguments

- `neighborhood`: a [`DispersalKernel`](@ref) based on any DynamicGrids.jl `Neighborhood`.

Pass grid name `Symbol`s to `R` and `W` type parameters to use specific grids.
"""
struct OutwardsDispersal{R,W,N<:AbstractKernel} <: SetNeighborhoodRule{R,W}
    neighborhood::N
end
function OutwardsDispersal{R,W}(; neighborhood=DispersalKernel{3}()) where {R,W} 
    OutwardsDispersal{R,W}(neighborhood)
end

@inline function applyrule!(data, rule::OutwardsDispersal{R,W}, state, I) where {R,W}
    state == zero(state) && return
    sum = zero(state)
    for (i, offset) in enumerate(offsets(rule))
        @inbounds propagules = state * kernel(rule)[i]
        @inbounds add!(data[W], propagules, I .+ offset...)
        sum += propagules
    end
    # Subtract from current cell, unless state is Bool
    state isa Bool || @inbounds sub!(data[W], sum, I...)
    return nothing
end

"""
    InwardsPopulationDispersal(neighborhood)
    InwardsPopulationDispersal(; neighborhood=DispersalKernel{3}())
    InwardsPopulationDispersal{R}(; neighborhood)
    InwardsPopulationDispersal{R,W}(; neighborhood)

Disperses to the current cells from the populations of the
surrounding cells, using a dispersal kernel deterministically.

This will make sense ecologically where cell populations are large, 
otherwise a randomised kernel may be more suitable.

The result should be identical to the matching [`OutwardsDispersal`](@ref).

# Arguments

- `neighborhood`: a [`DispersalKernel`](@ref) based on any DynamicGrids.jl `Neighborhood`.

Pass grid `Symbol`s to `R` or both `R` and `W` type parameters to use to specific grids.
"""
struct InwardsDispersal{R,W,N<:AbstractKernel} <: NeighborhoodRule{R,W}
    neighborhood::N
end
function InwardsDispersal{R,W}(; neighborhood=DispersalKernel{3}()) where {R,W}
    InwardsDispersal{R,W}(neighborhood)
end

@inline applyrule(data, rule::InwardsDispersal, N, I) = dot(neighborhood(rule))

"""
    InwardsPopulationDispersal(neighborhood)
    InwardsPopulationDispersal(; neighborhood=DispersalKernel{3}())
    InwardsPopulationDispersal{R}(; neighborhood)
    InwardsPopulationDispersal{R,W}(; neighborhood)

Implements deterministic dispersal from populations in neighboring cells to the current 
cell.

The result should be identical to those obtained substituting `InwardsDispersal` for 
[`OutwardsDispersal`](@ref) but will perform better when populations are spread across the
grid.


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

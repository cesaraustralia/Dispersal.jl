"""
    InwardsPopulationDispersal(neighborhood)
    InwardsPopulationDispersal(; neighborhood=DispersalKernel{3}())
    InwardsPopulationDispersal{R,W}(neighborhood)

Disperses to the current cells from the populations of the
surrounding cells, using a dispersal kernel deterministically.

This will only make sense where cell populations are large.

Pass grid name `Symbol`s to `R` and `W` type parameters to use specific grids.

The result should be identical to the matching [`OutwardsDispersal`](@ref).
"""
struct InwardsDispersal{R,W,N<:AbstractKernel} <: NeighborhoodRule{R,W}
    "Normalised proportions of dispersal to surrounding cells"
    neighborhood::N
end
InwardsDispersal{R,W}(; neighborhood=DispersalKernel{3}()) where {R,W} =
    InwardsDispersal{R,W}(neighborhood)

@inline applyrule(data, rule::InwardsDispersal, state, I) = 
    LinearAlgebra.dot(neighborhood(rule))

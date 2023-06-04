"""
    InwardsPopulationDispersal <: NeighborhoodRule

    InwardsPopulationDispersal(; kw...)
    InwardsPopulationDispersal{R}(; kw...)
    InwardsPopulationDispersal{R,W}(; kw...)

Implements deterministic dispersal from populations in neighboring cells to the current 
cell.

The result should be identical to those obtained substituting `InwardsDispersal` for 
[`OutwardsDispersal`](@ref) but will perform better when populations are spread across the
grid.

# Keywords

- `neighborhood`: Any DynamicGrids.jl `Neighborhood`, or an
    already constructed [`DispersalKernel`](@ref). Using this keyword means `radius` is
    ignored, and for a `DispersalKernel`, all other keywords are ignored.
- `neighborhood`: `Neighborhood` object specifying the range from the origin of the
    discretised dispersal kernal. Defaults to `Window(radius)`.
- `formulation`: kernel formulation object holding the exact form of the kernal.
    Default [`ExponentialKernel`](@ref).
- `cellsize`: the cell size of the discretised kernal (i.e. simulation grid size).
    Default is 1.0.
- `distancemethod`: [`DistanceMethod`](@ref) object for calculating distance between cells.
    The default is [`CentroidToCentroid`](@ref).

Pass grid `Symbol`s to `R` or both `R` and `W` type parameters to use to specific grids.
"""
struct InwardsDispersal{R,W,S<:Stencils.AbstractKernelStencil} <: NeighborhoodRule{R,W}
    stencil::S
end
function InwardsDispersal{R,W}(; kw...) where {R,W}
    InwardsDispersal{R,W}(DispersalKernel(; kw...))
end

@inline applyrule(data, rule::InwardsDispersal, N, I) = kernelproduct(rule)

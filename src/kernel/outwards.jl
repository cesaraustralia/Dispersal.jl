"""
    OutwardsPopulationDispersal <: SetNeighborhoodRule

    OutwardsPopulationDispersal(; kw...)
    OutwardsPopulationDispersal{R}(; kw...)
    OutwardsPopulationDispersal{R,W}(; kw...)

Implements deterministic dispersal from the current cell to populations in neighboring
cells.

This will make sense ecologically where cell populations are large,
otherwise a randomised kernel may be more suitable.

The result should be identical to those obtained substituting `OutwardsDispersal` for
[`InwardsDispersal`](@ref) but may be more efficient when a small number of cells are
occupied. Conversely, it will become less efficient when a large proportion of the grid
is occupied.

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

Pass grid name `Symbol`s to `R` and `W` type parameters to use specific grids.
"""
# Updated the OutwardsDispersal rule to include an optional mask
struct OutwardsDispersal{R,W,S<:Stencils.AbstractKernelStencil} <: SetNeighborhoodRule{R,W}
    stencil::S
    mask::M
end

function OutwardsDispersal{R,W}(stencil::S; mask=nothing) where {R,W,S<:Stencils.AbstractKernelStencil}
    OutwardsDispersal{R,W,S}(stencil, mask)
end

function OutwardsDispersal{R,W}(; kw...) where {R,W}
    OutwardsDispersal{R,W}(DispersalKernel(; kw...))
end

@inline function applyrule!(data, rule::OutwardsDispersal{R,W}, N, I) where {R,W}
    N == zero(N) && return nothing
    # Check if the current cell is masked, skip if it is
    if isnothing(rule.mask) || rule.mask[I...]
        return nothing
    end
    sum = zero(N)
    for (offset, k) in zip(offsets(rule), kernel(rule))
        target = I .+ offset
        (target_mod, inbounds) = _inbounds(Reflect(), size(data[1]), target...)
        if inbounds && (isnothing(rule.mask) || mask(data)[target_mod...])
            @inbounds propagules = N * k  
            @inbounds add!(data[W], propagules, target_mod...)  
            sum += propagules
        end
    end 
    @inbounds sub!(data[W], sum, I...)
    return nothing
end
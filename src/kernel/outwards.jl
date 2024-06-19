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
- `maskbehavior`: The default setting is `IgnoreMaskEdges()`. Use `CheckMaskEdges()` to indicate that the grid is 
    masked, enabling the rule to perform boundary checking at mask edges. Not using 
    `CheckMaskEdges()` on a masked grid may result in the loss of individuals at the edges, but it comes
    at a performance cost.

Pass grid name `Symbol`s to `R` and `W` type parameters to use specific grids.
"""

struct CheckMaskEdges end
struct IgnoreMaskEdges end

struct OutwardsDispersal{R,W,S<:Stencils.AbstractKernelStencil, M} <: SetNeighborhoodRule{R,W}
    stencil::S
    maskbehavior::M
end

# Constructors for OutwardsDispersal
function OutwardsDispersal{R,W}(stencil::S; maskbehavior::Union{CheckMaskEdges, IgnoreMaskEdges}=IgnoreMaskEdges()) where {R,W,S<:Stencils.AbstractKernelStencil}
    OutwardsDispersal{R,W,S,typeof(maskbehavior)}(stencil, maskbehavior)
end

function OutwardsDispersal{R,W}(; maskbehavior::Union{CheckMaskEdges, IgnoreMaskEdges}=IgnoreMaskEdges(), kw...) where {R,W}
    stencil = DispersalKernel(; kw...)
    OutwardsDispersal{R,W,typeof(stencil),typeof(maskbehavior)}(stencil, maskbehavior)
end

@inline function applyrule!(data, rule::OutwardsDispersal{R,W}, N, I) where {R,W}
    N == zero(N) && return nothing

    # Check if the current cell is masked, skip if it is
    mask_data = if rule.maskbehavior === IgnoreMaskEdges() nothing else DynamicGrids.mask(data) end
    if !isnothing(mask_data) && !mask_data[I...]
        return nothing
    end

    sum = zero(N)
    for (offset, k) in zip(offsets(rule), kernel(rule))
        target = I .+ offset
        (target_mod, inbounds) = DynamicGrids.inbounds(data, target)
        if inbounds && (isnothing(mask_data) || mask_data[target_mod...])
            @inbounds propagules = N * k  
            @inbounds add!(data[W], propagules, target_mod...)  
            sum += propagules
        end
    end
    @inbounds sub!(data[W], sum, I...)
    return nothing
end

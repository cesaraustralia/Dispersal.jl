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
- `mask_flag`: Use `Mask()` to apply masking or `NoMask()` to ignore masking. Default is `NoMask`.

Pass grid name `Symbol`s to `R` and `W` type parameters to use specific grids.
"""

struct Mask end
struct NoMask end

struct OutwardsDispersal{R,W,S<:Stencils.AbstractKernelStencil, M} <: SetNeighborhoodRule{R,W}
    stencil::S
    mask_flag::M
end

# Constructors for OutwardsDispersal
function OutwardsDispersal{R,W}(stencil::S; mask_flag::Union{Mask, NoMask}=NoMask()) where {R,W,S<:Stencils.AbstractKernelStencil}
    OutwardsDispersal{R,W,S,typeof(mask_flag)}(stencil, mask_flag)
end

function OutwardsDispersal{R,W}(; mask_flag::Union{Mask, NoMask}=NoMask(), kw...) where {R,W}
    stencil = DispersalKernel(; kw...)
    OutwardsDispersal{R,W,typeof(stencil),typeof(mask_flag)}(stencil, mask_flag)
end

# @inline function applyrule!(data, rule::OutwardsDispersal{R,W}, N, I) where {R,W}
#     N == zero(N) && return nothing
#     mask_data = rule.mask_flag === NoMask() ? nothing : DynamicGrids.mask(data)
#     sum = zero(N)

#     if isnothing(mask_data)
#         # If there is no mask
#         for (offset, k) in zip(offsets(rule), kernel(rule))
#             @inbounds propagules = N * k
#             @inbounds add!(data[W], propagules, I .+ offset...)
#             sum += propagules
#         end
#     elseif !mask_data[I...]
#         # If there is a mask and the source cell is masked
#         return nothing
#     else
#         for (offset, k) in zip(offsets(rule), kernel(rule))
#             (target_mod, inbounds) = inbounds(data, I .+ offset)
#             if inbounds && mask_data[target_mod...]
#                 @inbounds propagules = N * k  
#                 @inbounds add!(data[W], propagules, target_mod...)  
#                 sum += propagules
#             end
#         end
#     end

#     @inbounds sub!(data[W], sum, I...)
#     return nothing
# end

@inline function applyrule!(data, rule::OutwardsDispersal{R,W}, N, I) where {R,W}
    N == zero(N) && return nothing

    # Check if the current cell is masked, skip if it is
    mask_data = if rule.mask_flag === NoMask() nothing else DynamicGrids.mask(data) end
    if !isnothing(mask_data) && !mask_data[I...]
        return nothing
    end
    
    sum = zero(N)
    for (offset, k) in zip(offsets(rule), kernel(rule))
        target = I .+ offset
        inbounds = if isnothing(mask_data)    
            true
        else        
            (target_mod, inbounds) = DynamicGrids.inbounds(data, target)
            mask_data[target_mod...]
        end
        if inbounds
            propagules = N * k  
            @inbounds add!(data[W], propagules, target...)  
            sum += propagules
        end
    end
    @inbounds sub!(data[W], sum, I...)
    return nothing
end
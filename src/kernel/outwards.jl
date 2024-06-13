"""
    OutwardsPopulationDispersal <: SetNeighborhoodRule

    OutwardsPopulationDispersal(; kw...)
    OutwardsPopulationDispersal{R}(; kw...)
    OutwardsPopulationDispersal{R,W}(; kw...)

Implements deterministic dispersal from the current cell to populations in neighboring
cells. ups

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
struct _Mask end
const NoMask = nothing

struct OutwardsDispersal{R,W,S<:Stencils.AbstractKernelStencil, M} <: SetNeighborhoodRule{R,W}
    stencil::S
    mask_flag::M
end

# Constructors for OutwardsDispersal
function OutwardsDispersal{R,W}(stencil::S; mask_flag=NoMask) where {R,W,S<:Stencils.AbstractKernelStencil}
    OutwardsDispersal{R,W,S,typeof(mask_flag)}(stencil, mask_flag)
end

function OutwardsDispersal{R,W}(; mask_flag=NoMask, kw...) where {R,W}
    stencil = DispersalKernel(; kw...)
    OutwardsDispersal{R,W,typeof(stencil),typeof(mask_flag)}(stencil, mask_flag)
end

@inline function applyrule!(data, rule::OutwardsDispersal{R,W}, N, I) where {R,W}
    N == zero(N) && return nothing
    
    # Check if the current cell is masked, skip if it is
    mask_data = if rule.mask_flag === NoMask nothing else mask(data) end
    if !isnothing(mask_data) && !mask_data[I...]
        return nothing
    end
    
    sum = zero(N)
    for (offset, k) in zip(offsets(rule), kernel(rule))
        target = I .+ offset
        (target_mod, inbounds) = inbounds(data, target)
        if inbounds && (isnothing(mask_data) || mask_data[target_mod...])
            @inbounds propagules = N * k  
            @inbounds add!(data[W], propagules, target_mod...)  
            sum += propagules
        end
    end
    @inbounds sub!(data[W], sum, I...)
    return nothing
end

# @inline function applyrule(data, rule::OutwardsDispersal{R,W}, N, I) where {R,W}
#     applyrule!(data, rule, N, I)
# end

################# TESTING #################
# Define a mask
mask_data = fill(true, 10, 10)
for i in 3:9 
    mask_data[2, i] = false
    mask_data[3, i] = false
    mask_data[4, i] = false
    mask_data[5, i] = false
    mask_data[6, i] = false
    mask_data[7, i] = false
    mask_data[8, i] = false
end

# Create OutwardsDispersal with and without a mask
outdisp_with_mask = OutwardsDispersal(
    formulation=ExponentialKernel(λ=0.0125),
    distancemethod=AreaToArea(30),
    mask_flag=_Mask()
)

outdisp_without_mask = OutwardsDispersal(
    formulation=ExponentialKernel(λ=0.0125),
    distancemethod=AreaToArea(30)
)

# Initialize the simulation grid
init = fill(100.0, (10,10))
init[5,5] = 0.0

# Create ruleset and outputs
rule_with_mask = Ruleset(outdisp_with_mask; boundary=Reflect())
rule_without_mask = Ruleset(outdisp_without_mask; boundary=Reflect())

# Run the simulation with a mask
out_with_mask = ArrayOutput(init; tspan=1:100, mask=mask_data)
a = sim!(out_with_mask, rule_with_mask)
a[1][5,5]
a[2][5,5]
# Run the simulation without a mask
out_without_mask = ArrayOutput(init; tspan=1:100)
sim!(out_without_mask, rule_without_mask)


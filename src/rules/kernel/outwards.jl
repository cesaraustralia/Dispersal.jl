
"""
Abstract supertype that extends `ManualNeighborhoodRule` for neighborhood-based
dispersal rules that update surounding cells based on the values of the
current cell, as if dispersing outwards from the current cell.

The result should be identical to [`InwardsDispersal`](@ref) but may be more
efficient than when a small number of cells are occupied. It is less efficient
when a large proportion of the grid is occupied.
"""
abstract type OutwardsDispersal{R,W} <: ManualNeighborhoodRule{R,W} end

# The same rule is used for Binary and Population rules,
# with different setneighbour! methods

@inline function applyrule!(data, rule::OutwardsDispersal{R,W}, state, cellindex) where {R,W}
    state == zero(state) && return
    hood = neighborhood(rule)
    sum = mapsetneighbor!(data[W], hood, rule, state, cellindex)
    # Subtract from current cell, unless state is Bool
    state isa Bool || @inbounds return add!(data[W], -sum, cellindex...)
    return nothing
end

"""
    OutwardsBinaryDispersal(neighborhood)
    OutwardsBinaryDispersal(; neighborhood=DispersalKernel{3}())
    OutwardsBinaryDispersal{R,W}(neighborhood)
    OutwardsBinaryDispersal{R,W}(; neighborhood=DispersalKernel{3}())

Cells in the surrounding neighborhood have some propability of
invasion if the current cell is occupied.

Pass grid name `Symbol`s to `R` and `W` type parameters to use specific grids.
"""
struct OutwardsBinaryDispersal{R,W,NH,PT} <: OutwardsDispersal{R,W}
    "Normalised proportions of dispersal to surrounding cells"
    neighborhood::NH
    "A real number between one and zero"
    prob_threshold::PT
end
OutwardsBinaryDispersal{R,W}(;
    neighborhood=DispersalKernel{3}(),
    prob_threshold=Param(0.1, bounds=(0.0, 1.0)),
) where {R,W} = OutwardsBinaryDispersal{R,W}(neighborhood, prob_threshold)

@inline function setneighbor!(
    data::WritableGridData, hood::Neighborhood, rule::OutwardsBinaryDispersal,
    state::Bool, hood_index, dest_index
)
    @inbounds rand() * kernel(hood)[hood_index...] > rule.prob_threshold || return zero(state)
    @inbounds or!(data, oneunit(state), dest_index...)
    return oneunit(state)
end

"""
    OutwardsPopulationDispersal(neighborhood)
    OutwardsPopulationDispersal(; neighborhood=DispersalKernel{3}())
    OutwardsPopulationDispersal{R,W}(neighborhood)
    OutwardsPopulationDispersal{R,W}(; neighborhood=DispersalKernel{3}())

Dispersal reduces the current cell population, increasing the populations of the
cells in the surrounding neighborhood.

Pass grid name `Symbol`s to `R` and `W` type parameters to use specific grids.
"""
struct OutwardsPopulationDispersal{R,W,NH} <: OutwardsDispersal{R,W}
    "Normalised proportions of dispersal to surrounding cells"
    neighborhood::NH
end
OutwardsPopulationDispersal{R,W}(;
    neighborhood=DispersalKernel{3}(),
) where {R,W} = OutwardsPopulationDispersal{R,W}(neighborhood)


@inline function setneighbor!(
    data::WritableGridData, hood::Neighborhood, rule::OutwardsPopulationDispersal,
    state::AbstractFloat, hood_index, dest_index
)
    @inbounds propagules = state * kernel(hood)[hood_index...]
    @inbounds add!(data, propagules, dest_index...)
    return propagules
end

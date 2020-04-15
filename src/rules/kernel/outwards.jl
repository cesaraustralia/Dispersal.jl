
"""
Abstract supertype that extends `PartialNeighborhoodRule` for neighborhood-based 
dispersal rules that update surounding cells based on the values of the 
current cell, as if dispersing outwards from the current cell.

The result should be identical to [`InwardsDispersal`](@ref) but may be more 
efficient than when a small number of cells are occupied. It is less efficient 
when a large proportion of the grid is occupied.
"""
abstract type OutwardsDispersal{R,W} <: PartialNeighborhoodRule{R,W} end

"""
    OutwardsBinaryDispersal(neighborhood)
    OutwardsBinaryDispersal(; neighborhood=DispersalKernel{3}())
    OutwardsBinaryDispersal{R,W}(neighborhood)

Cells in the surrounding neighborhood have some propability of 
invasion if the current cell is occupied.

Pass grid name `Symbol`s to `R` and `W` type parameters to use specific grids.

$(FIELDDOCTABLE)
"""
@Kernel @Probabilistic struct OutwardsBinaryDispersal{R,W} <: OutwardsDispersal{R,W} end

"""
    OutwardsPopulationDispersal(neighborhood)
    OutwardsPopulationDispersal(; neighborhood=DispersalKernel{3}())
    OutwardsPopulationDispersal{R,W}(neighborhood)

Dispersal reduces the current cell population, increasing the populations of the 
cells in the surrounding neighborhood.

Pass grid name `Symbol`s to `R` and `W` type parameters to use specific grids.

$(FIELDDOCTABLE)
"""
@Kernel struct OutwardsPopulationDispersal{R,W} <: OutwardsDispersal{R,W} end

@inline applyrule!(rule::OutwardsDispersal{R,W}, data, state, index) where {R,W} = begin
    state == zero(state) && return
    hood = neighborhood(rule)
    sum = mapsetneighbor!(data[W], hood, rule, state, index)
    update_state!(data[W], hood, state, index, sum)
    return
end

@inline update_state!(grid, hood, state::AbstractFloat, index, sum) = 
    grid[index...] = state - sum
@inline update_state!(grid, hood, state, index, sum) = state

@inline setneighbor!(data, hood, rule::OutwardsPopulationDispersal, 
                     state::AbstractFloat, hood_index, dest_index) = begin
    @inbounds propagules = state * kernel(hood)[hood_index...]
    @inbounds data[dest_index...] += propagules
    propagules
end

@inline setneighbor!(data, hood, rule::OutwardsBinaryDispersal,
                     state::Integer, hood_index, dest_index) = begin
    @inbounds rand() * kernel(hood)[hood_index...] > rule.prob_threshold || return zero(state)
    @inbounds data[dest_index...] += oneunit(state)
    oneunit(state)
end

@inline setneighbor!(data, hood, rule::OutwardsBinaryDispersal, 
                     state::Bool, hood_index, dest_index) = begin
    @inbounds rand() * kernel(hood)[hood_index...] > rule.prob_threshold || return zero(state)
    @inbounds data[dest_index...] |= oneunit(state)
    oneunit(state)
end

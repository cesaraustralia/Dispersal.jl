
"""
Extends PartialNeighborhoodRule for outwards dispersal.

Outwards dispersal calculates dispersal *from* the current cell to cells
in its [`DispersalNeighborhood`](@ref). This should be more efficient than 
inwards dispersal when a small number of cells are occupied, but less efficient 
when a large proportion of the grid is occupied.
"""
abstract type OutwardsDispersal{R} <: PartialNeighborhoodRule{R} end

# Get the radius from the kernel for all OutwardsDispersal
(::Type{T})(kernel, args...) where T <: OutwardsDispersal = 
    T{radius(kernel),typeof(kernel),typeof.(args)...}(kernel, args...)

"""
Cells in the surrounding [`DispersalNeighborhood`](@ref) have some propability of 
invasion if the current cell is occupied.
$(FIELDDOCTABLE)
"""
@Kernel @Probabilistic struct OutwardsBinaryDispersal{R} <: OutwardsDispersal{R} end

"""
Dispersal reduces the current cell population, increasing the populations of the 
cells in the surrounding [`DispersalNeighborhood`](@ref).
$(FIELDDOCTABLE)
"""
@Kernel struct OutwardsPopulationDispersal{R} <: OutwardsDispersal{R} end

DynamicGrids.radius(rule::OutwardsDispersal) = radius(rule.neighborhood)


@inline applyrule!(rule::OutwardsDispersal, data, state, index) = begin
    hood = neighborhood(rule)
    sum = mapreduceneighbors(setneighbor!, data, hood, rule, state, index)
    update_state!(data, hood, state, index, sum)
end

@inline update_state!(data, hood, state::AbstractFloat, index, sum) = 
    data[index...] = state - sum
@inline update_state!(data, hood, state, index, sum) = state

@inline setneighbor!(data, hood, rule::OutwardsPopulationDispersal, state::AbstractFloat, hood_index, dest_index) = begin
    @inbounds propagules = state * kernel(hood)[hood_index...]
    @inbounds data[dest_index...] += propagules
    propagules
end

@inline setneighbor!(data, hood, rule::OutwardsBinaryDispersal, state::Integer, hood_index, dest_index) = begin
    @inbounds rand() * kernel(hood)[hood_index...] > rule.prob_threshold || return zero(state)
    @inbounds data[dest_index...] += oneunit(state)
    oneunit(state)
end

@inline setneighbor!(data, hood, rule::OutwardsBinaryDispersal, state::Bool, hood_index, dest_index) = begin
    @inbounds rand() * kernel(hood)[hood_index...] > rule.prob_threshold || return zero(state)

    @inbounds data[dest_index...] |= oneunit(state)
    oneunit(state)
end

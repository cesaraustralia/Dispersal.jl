
"""
Extends AbstractPartialNeighborhoodRule for outwards dispersal.

Outwards dispersal calculates dispersal *from* the current cell to cells
in its [`DispersalNeighborhood`](@ref). This should be more efficient than 
inwards dispersal when a small number of cells are occupied, but less efficient 
when a large proportion of the grid is occupied.
"""
abstract type AbstractOutwardsDispersal{R} <: AbstractPartialNeighborhoodRule{R} end

# Get the radius from the kernel for all AbstractOutwardsDispersal
(::Type{T})(kernel, args...) where T <: AbstractOutwardsDispersal = 
    T{radius(kernel),typeof(kernel),typeof.(args)...}(kernel, args...)

"""
Cells in the surrounding [`DispersalNeighborhood`](@ref) have some propability of 
invasion if the current cell is occupied.
$(FIELDDOCTABLE)
"""
@Kernel @Probabilistic struct OutwardsBinaryDispersal{R} <: AbstractOutwardsDispersal{R} end

"""
Dispersal reduces the current cell population, increasing the populations of the 
cells in the surrounding [`DispersalNeighborhood`](@ref).
$(FIELDDOCTABLE)
"""
@Kernel struct OutwardsPopulationDispersal{R} <: AbstractOutwardsDispersal{R} end

DynamicGrids.radius(rule::AbstractOutwardsDispersal) = radius(rule.neighborhood)


@inline applyrule!(rule::AbstractOutwardsDispersal, data, state, index) = begin
    hood = rule.neighborhood
    sum = mapreduceneighbors(setneighbor!, data, hood, rule, state, index)
    update_state!(data, hood, state, index, sum)
end

@inline update_state!(data, hood, state::AbstractFloat, index, sum) = 
    data[index...] = state - sum
@inline update_state!(data, hood, state, index, sum) = state

@inline setneighbor!(data, hood, rule::OutwardsPopulationDispersal, state::AbstractFloat, hood_index, dest_index) = begin
    @inbounds propagules = state * hood.kernel[hood_index...]
    @inbounds data[dest_index...] += propagules
    propagules
end

@inline setneighbor!(data, hood, rule::OutwardsBinaryDispersal, state::Integer, hood_index, dest_index) = begin
    @inbounds rand() * hood.kernel[hood_index...] > rule.prob_threshold || return zero(state)
    @inbounds data[dest_index...] += oneunit(state)
    oneunit(state)
end

@inline setneighbor!(data, hood, rule::OutwardsBinaryDispersal, state::Bool, hood_index, dest_index) = begin
    @inbounds rand() * hood.kernel[hood_index...] > rule.prob_threshold || return zero(state)

    @inbounds data[dest_index...] |= oneunit(state)
    oneunit(state)
end

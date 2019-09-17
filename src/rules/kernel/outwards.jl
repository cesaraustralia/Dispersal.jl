
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
    neighbors(rule.neighborhood, rule, data, state, index)
    data[index...]
end

@inline neighbors(hood, rule::AbstractPartialNeighborhoodRule, data, state, index) = begin
    r = radius(hood)
    propagules = zero(state)
    # Loop over dispersal kernel grid dimensions
    for x = one(r):2r + one(r)
        xs = x + index[2] - r - one(r)
        @simd for y = one(r):2r + one(r)
            ys = y + index[1] - r - one(r)
            # Update cumulative value, and cell value for outwards dispersal
            propagules += update_cell!(hood, rule, data, state, (y, x), (ys, xs))
        end
    end
    update_state(rule, data, state, index, propagules)
    propagules
end

update_state(rule, data, state::AbstractFloat, index, propagules) = 
    data[index...] -= propagules
update_state(rule, data, state, index, propagules) = nothing

@inline update_cell!(hood, rule, data, state::AbstractFloat, hood_index, dest_index) = begin
    @inbounds propagules = state * hood.kernel[hood_index...]
    @inbounds data[dest_index...] += propagules
    propagules
end

@inline update_cell!(hood, rule, data, state::Integer, hood_index, dest_index) = begin
    @inbounds rand() * hood.kernel[hood_index...] > rule.prob_threshold || return zero(state)

    @inbounds data[dest_index...] += oneunit(state)
    oneunit(state)
end

@inline update_cell!(hood, rule, data, state::Bool, hood_index, dest_index) = begin
    @inbounds rand() * hood.kernel[hood_index...] > rule.prob_threshold || return zero(state)

    @inbounds data[dest_index...] |= oneunit(state)
    oneunit(state)
end

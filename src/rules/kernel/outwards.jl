
"""
Extends AbstractPartialNeighborhoodRule for outwards dispersal
"""
abstract type AbstractOutwardsDispersal <: AbstractPartialNeighborhoodRule end

"""
Binary binary dispersal within a [`DispersalNeighborhood`](@ref)

Outwards dispersal calculates dispersal *from* the current cell to cells
in its neighborhood. This should be more efficient than inwards
dispersal when a small number of cells are occupied, but less efficient when a large
proportion of the grid is occupied.

Surrounding cells are invaded if the current cell is occupied.
$(FIELDDOCTABLE)
"""
@Probabilistic @Kernel struct OutwardsBinaryDispersal{} <: AbstractOutwardsDispersal end

"""
Disperses from the current cells population to the populations of the surrounding cells,
using a dispersal kernel.
$(FIELDDOCTABLE)
"""
@Kernel struct OutwardsPopulationDispersal{} <: AbstractOutwardsDispersal end

CellularAutomataBase.radius(rule::AbstractOutwardsDispersal) = radius(rule.neighborhood)


@inline applyrule!(rule::AbstractOutwardsDispersal, data, state, index, args...) = begin
    state == zero(state) && return state # Ignore empty cells
    neighbors(rule.neighborhood, rule, data, state, index, args...)
    data.dest[index...]
end

@inline neighbors(hood, rule::AbstractPartialNeighborhoodRule, data, state, index, args...) = begin
    r = radius(hood)
    propagules = zero(state)
    # Loop over dispersal kernel grid dimensions
    for x = one(r):2r + one(r)
        xs = x + index[2] - r - one(r)
        @simd for y = one(r):2r + one(r)
            ys = y + index[1] - r - one(r)
            # Update cumulative value, and cell value for outwards dispersal
            propagules += update_cell!(hood, rule, data, state, (y, x), (ys, xs), args...)
        end
    end
    update_state(rule, data, state, index, propagules)
    propagules
end

update_state(rule, data, state::AbstractFloat, index, propagules) = 
    data.dest[index...] -= propagules
update_state(rule, data, state, index, propagules) = nothing

@inline update_cell!(hood, rule, data, state::AbstractFloat, hood_index, dest_index, args...) = begin
    @inbounds propagules = state * hood.kernel[hood_index...]
    @inbounds data.dest[dest_index...] += propagules
    propagules
end

@inline update_cell!(hood, rule, data, state::Integer, hood_index, dest_index, args...) = begin
    @inbounds rand() * hood.kernel[hood_index...] > rule.prob_threshold || return zero(state)

    @inbounds data.dest[dest_index...] += oneunit(state)
    oneunit(state)
end

@inline update_cell!(hood, rule, data, state::Bool, hood_index, dest_index, args...) = begin
    @inbounds rand() * hood.kernel[hood_index...] > rule.prob_threshold || return zero(state)

    @inbounds data.dest[dest_index...] |= oneunit(state)
    oneunit(state)
end

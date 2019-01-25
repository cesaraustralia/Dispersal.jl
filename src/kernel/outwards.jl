
"Extend to modify [`OutwardsBinaryDispersal`](@ref)"
abstract type AbstractOutwardsDispersal <: AbstractPartialNeighborhoodModel end

"""
Binary binary dispersal within a [`DispersalNeighborhood`](@ref)

Outwards dispersal calculates dispersal *from* the current cell to cells
in its neighborhood. This should be more efficient than inwards
dispersal when a small number of cells are occupied, but less efficient when a large
proportion of the grid is occupied.
"""
@Probabilistic @Kernel struct OutwardsBinaryDispersal{} <: AbstractOutwardsDispersal end

@Kernel struct OutwardsPopulationDispersal{} <: AbstractOutwardsDispersal end

Cellular.radius(model::AbstractOutwardsDispersal) = radius(model.neighborhood)
Cellular.temp_neighborhood(model::AbstractOutwardsDispersal) = temp_neighborhood(model.neighborhood)

"""
    rule(model::AbstractOutwardsDispersal, data, state, index, args...)
Runs rule for of [`AbstractOutwardsDispersal`](@ref) dispersal.

Surrounding cells are invaded if the current cell is occupied and they have
suitable habitat. Otherwise they keeps their current state.
"""
@inline rule!(model::AbstractOutwardsDispersal, data, state, index, args...) = begin
    state == zero(state) && return state # Ignore empty cells
    neighbors(model.neighborhood, model, data, state, index, args...)
    data.dest[index...]
end

@inline neighbors(hood, model::AbstractPartialNeighborhoodModel, data, state, index, args...) = begin
    r = hood.radius
    propagules = zero(state)
    # Loop over dispersal kernel grid dimensions
    for x = one(r):2r + one(r)
        xs = x + index[2] - r - one(r)
        @simd for y = one(r):2r + one(r)
            ys = y + index[1] - r - one(r)
            # Update cumulative value, and cell value for outwards dispersal
            propagules += update_cell!(hood, model, data, state, (y, x), (ys, xs), args...)
        end
    end
    update_state(model, data, state, index, propagules)
    propagules
end

update_state(model, data, state::AbstractFloat, index, propagules) = data.dest[index...] -= propagules
update_state(model, data, state, index, propagules) = nothing

@inline update_cell!(hood, model, data, state::AbstractFloat,
                   hood_index, dest_index, args...) = begin
    @inbounds propagules = state * hood.kernel[hood_index...]
    @inbounds data.dest[dest_index...] += propagules
    propagules
end

@inline update_cell!(hood, model, data, state::Integer,
                   hood_index, dest_index, args...) = begin
    @inbounds rand() * hood.kernel[hood_index...] > model.prob_threshold || return zero(state)

    @inbounds data.dest[dest_index...] += oneunit(state)
    oneunit(state)
end

@inline update_cell!(hood, model, data, state::Bool,
                   hood_index, dest_index, args...) = begin
    @inbounds rand() * hood.kernel[hood_index...] > model.prob_threshold || return zero(state)

    @inbounds data.dest[dest_index...] |= oneunit(state)
    oneunit(state)
end


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

@Fraction @Kernel struct OutwardsPopulationDispersal{} <: AbstractOutwardsDispersal end


"""
    rule(model::AbstractOutwardsDispersal, data, state, index, args...)
Runs rule for of [`AbstractOutwardsDispersal`](@ref) dispersal.

Surrounding cells are invaded if the current cell is occupied and they have
suitable habitat. Otherwise they keeps their current state.
"""
@inline rule!(model::AbstractOutwardsDispersal, data, state, index, args...) = begin
    state == zero(state) && return state # Ignore empty cells
    propagules = neighbors(model.neighborhood, model, data, state, index, args...)
    data.dest[index...] -= propagules
end

@inline neighbors(hood, model::AbstractPartialNeighborhoodModel, data, state, index, args...) = begin
    r = hood.radius
    cc = zero(state)
    # Loop over dispersal kernel grid dimensions
    for x = one(r):2r + one(r)
        xs = x + index[2] - r - one(r)
        @simd for y = one(r):2r + one(r)
            ys = y + index[1] - r - one(r)
            # Update cumulative value, and cell value for outwards dispersal
            cc += update_cell!(hood, model, data, state, (y, x), (ys, xs), args...)
        end
    end
    cc
end

@inline update_cell!(hood, model, data, state::AbstractFloat,
                   hood_index, dest_index, args...) = begin
    @inbounds propagules = state * hood.kernel[hood_index...]
    @inbounds data.dest[dest_index...] += propagules
    propagules
end

@inline update_cell!(hood, model, data, state::Integer,
                   hood_index, dest_index, args...) = begin
    @inbounds spec_rand(source, Float64, args...) * hood.kernel[hood_index...] > model.prob_threshold || return zero(state)

    @inbounds data.dest[dest_index...] += oneunit(state)
    oneunit(state)
end

@inline update_cell!(hood, model, data, state::Bool,
                   hood_index, dest_index, args...) = begin
    @inbounds spec_rand(source, Float64, args...) * hood.kernel[hood_index...] > model.prob_threshold || return zero(state)

    @inbounds data.dest[dest_index...] |= oneunit(state)
    oneunit(state)
end

using LinearAlgebra

import Cellular: radius

Cellular.radius(hood::DispersalNeighborhood) = hood.radius

"""
    neighbors(hood::DispersalNeighborhood, state, index, t, source, dest, args...)

Returns nieghbors for a [`DispersalNeighborhood`](@ref), looping over
the array of dispersal propabilities.
"""
@inline neighbors(hood, model::AbstractNeighborhoodModel, data, state, index, args...) =
    @inbounds return dot(data.modelmem, hood.kernel)

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
                   hood_index, dest_index, layers, args...) = begin
    @inbounds propagules = state * model.fraction * hood.kernel[hood_index...]
    @inbounds data.dest[dest_index...] += propagules
    propagules
end

@inline update_cell!(hood, model, data, state::Integer,
                   hood_index, dest_index, layers, args...) = begin
    @inbounds spec_rand(source, Float64, args...) * hood.kernel[hood_index...] > model.prob_threshold || return zero(state)

    @inbounds data.dest[dest_index...] += oneunit(state)
    oneunit(state)
end

@inline update_cell!(hood, model, data, state::Bool,
                   hood_index, dest_index, layers, args...) = begin
    @inbounds spec_rand(source, Float64, args...) * hood.kernel[hood_index...] > model.prob_threshold || return zero(state)

    @inbounds data.dest[dest_index...] |= oneunit(state)
    oneunit(state)
end

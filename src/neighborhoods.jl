"""
    neighbors(hood::DispersalNeighborhood, state, index, t, source, dest, args...)

Returns nieghbors for a [`DispersalNeighborhood`](@ref), looping over
the array of dispersal propabilities.
"""
neighbors(hood::DispersalNeighborhood, model, data, state, index, args...) = begin
    cc = zero(state)
    r = hood.radius

    # Loop over dispersal kernel grid dimensions
    for b = -r:r, a = -r:r

        # Check boundaries
        source_index, is_inbounds = inbounds(index .+ (b, a), size(data.source), hood.overflow)
        is_inbounds || continue
        hood_index = (b, a) .+ (r + one(r))

        # Update cumulative value, and cell value for outwards dispersal
        cc += update_cell!(hood, model, data, state, hood_index, source_index, args...)
    end
    return cc
end

update_cell!(hood::DispersalNeighborhood{:outwards}, model, data, state::AbstractFloat,
                   hood_index, dest_index, layers, args...) = begin
    @inbounds propagules = state * model.fraction * hood.kernel[hood_index...]
    @inbounds data.dest[dest_index...] += propagules
    propagules
end

update_cell!(hood::DispersalNeighborhood{:outwards}, model, data, state::Integer,
                   hood_index, dest_index, layers, args...) = begin
    spec_rand(source, Float64, args...) * hood.kernel[hood_index...] > model.prob_threshold || return zero(state)

    @inbounds data.dest[dest_index...] = oneunit(state)
    oneunit(state)
end

update_cell!(hood::DispersalNeighborhood{:inwards}, model, data, state, hood_index, 
             source_index, args...) = begin
    @inbounds data.source[source_index...] * hood.kernel[hood_index...]
end

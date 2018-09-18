""" neighbors(hood::DispersalNeighborhood, state, row, col, t, source, dest, args...)
Returns nieghbors for a [`DispersalNeighborhood`](@ref), looping over
the array of dispersal propabilities.
"""
neighbors(hood::DispersalNeighborhood, model, state, row, col, t, source, dest, args...) = begin
    cc = zero(state) 
    r = hood.radius

    # Loop over dispersal kernel grid dimensions
    # TODO use an OffsetArray with 0,0 as the center
    for b = -r:r, a = -r:r
        # Check boundaries
        y, x, is_inbounds = inbounds((row, col) .+ (b, a), size(source), hood.overflow)
        is_inbounds || continue
        # Update cumulative value, and cell value for outwards dispersal
        cc += update_cell!(hood, model, state, t, source, dest, (b, a) .+ (r + one(r)), (y, x), args...)
    end
    return cc
end

update_cell!(hood::DispersalNeighborhood{:outwards}, model, state::AbstractFloat,
                   t, source, dest, hood_index, dest_index, layers, args...) = begin
    # Invade the cell
    # TODO: randomisation logic
    @inbounds dest[dest_index...] = spec_rand(source, Float64, args...) * hood.kernel[hood_index...]
end

update_cell!(hood::DispersalNeighborhood{:outwards}, model, state::Integer,
                   t, source, dest, hood_index, dest_index, layers, args...) = begin
    suitability(layers, dest_index, t) > model.suitability_threshold || return zero(state)
    spec_rand(source, Float64, args...) * hood.kernel[hood_index...] > model.prob_threshold || return zero(state)
    # Invade the cell
    @inbounds dest[dest_index...] = oneunit(state)
end

update_cell!(hood::DispersalNeighborhood{:inwards}, model, state, t, source, dest, 
             hood_index, source_index, args...) = begin
    source[source_index...] * hood.kernel[hood_index...]
end

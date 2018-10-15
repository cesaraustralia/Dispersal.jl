"""
    neighbors(hood::DispersalNeighborhood, state, row, col, t, source, dest, args...)

Returns nieghbors for a [`DispersalNeighborhood`](@ref), looping over
the array of dispersal propabilities.
"""
neighbors(hood::DispersalNeighborhood, model, state, row, col, t, source, dest, args...) = begin
    cc = zero(state)
    r = hood.radius

    # Loop over dispersal kernel grid dimensions
    for b = -r:r, a = -r:r

        # Check boundaries
        source_index, is_inbounds = inbounds((row, col) .+ (b, a), size(source), hood.overflow)
        is_inbounds || continue
        hood_index = (b, a) .+ (r + one(r))

        # Update cumulative value, and cell value for outwards dispersal
        cc += update_cell!(hood, model, state, t, source, dest, hood_index, source_index, args...)
    end
    return cc
end

update_cell!(hood::DispersalNeighborhood{:outwards}, model, state::AbstractFloat,
                   t, source, dest, hood_index, dest_index, layers, args...) = begin
    @inbounds propagules = state * model.fraction * hood.kernel[hood_index...]
    @inbounds dest[dest_index...] += propagules
    propagules
end

update_cell!(hood::DispersalNeighborhood{:outwards}, model, state::Integer,
                   t, source, dest, hood_index, dest_index, layers, args...) = begin
    spec_rand(source, Float64, args...) * hood.kernel[hood_index...] > model.prob_threshold || return zero(state)

    @inbounds dest[dest_index...] = oneunit(state)
    oneunit(state)
end

update_cell!(hood::DispersalNeighborhood{:inwards}, model, state, t, source, dest,
             hood_index, source_index, args...) = begin
    @inbounds source[source_index...] * hood.kernel[hood_index...]
end

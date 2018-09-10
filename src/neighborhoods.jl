"""
    DispersalNeighborhood(; f=d -> exponential(d, 1), radius=3, overflow=Skip())
Constructor for neighborhoods, using a dispersal kernel function and a cell radius.

### Keyword Arguments:
- `f::Function`: any function that accepts a Number argument and returns a propbability between 0.0 and 1.0
- `radius::Integer`: a positive integer
- `overflow = Skip()
"""
build_dispersal_kernel(f, params, init, cellsize, r) = begin
    params = typeof(params) <: Tuple ? params : (params,)
    sze = 2r + one(r)
    kernel = zeros(Float64, sze, sze)
    # Paper: l. 97
    for y = -r:r, x = -r:r
        kernel[y+r+one(y), x+r+one(x)] = f(sqrt(y^2 + x^2) * cellsize, params...)
    end
    # Normalise
    kernel ./= sum(kernel)
    SMatrix{sze,sze}(kernel)
end

# Paper: l. 96
exponential(d, a) = exp(-d / a)

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
                   t, source, dest, hood_index, dest_index, args...) = begin
    # Invade the cell
    # TODO: randomisation logic
    dest[dest_index...] = rand() * hood.kernel[hood_index...]
end

update_cell!(hood::DispersalNeighborhood{:outwards}, model, state::Integer,
                   t, source, dest, hood_index, dest_index, layers, args...) = begin
    suitability(layers, dest_index, t) > model.suitability_threshold || return zero(state)
    rand() * hood.kernel[hood_index...] > model.prob_threshold || return zero(state)
    # Invade the cell
    dest[dest_index...] = oneunit(state)
end

update_cell!(hood::DispersalNeighborhood{:inwards}, model, state, t, source, dest, 
             hood_index, source_index, args...) = begin
    @inbounds source[source_index...] * hood.kernel[hood_index...] 
end

using LinearAlgebra

"""
    neighbors(hood::DispersalNeighborhood, state, index, t, source, dest, args...)

Returns nieghbors for a [`DispersalNeighborhood`](@ref), looping over
the array of dispersal propabilities.
"""
neighbors(hood::DispersalNeighborhood, model, data, state, index, args...) = begin
    r = hood.radius

    y1 = -r
    for n = -r:0
        n + index[1] > 0 && break
        y1 = n
    end
    s = size(data.source, 1)
    y2 = r
    for n = r:1
        n + index[1] < s && break
        y2 = n
    end

    x1 = -r
    for n = -r:0
        n + index[2] > 0 && break
        x1 = n
    end
    s = size(data.source, 2)
    x2 = r
    for n = r:1
        n + index[2] < s && break
        x2 = n
    end

    cc = zero(state)
    # Loop over dispersal kernel grid dimensions
    for x = x1:x2
        xh = x + r + one(r)
        xs = x + index[2]
        @inbounds @simd for y = y1:y2
            yh = y + r + one(r)
            ys = y + index[1]
            # Update cumulative value, and cell value for outwards dispersal
            cc += update_cell!(hood, model, data, state, (yh, xh), (ys, xs), args...)
        end
    end
    cc
end


@inline neighbors(hood::DispersalNeighborhood{:inwards}, model, data, state, index, args...) = begin
    @inbounds dot(data.modeldata[1].loc, hood.kernel)
    # r = hood.radius
    # sum = zero(eltype(hood.kernel))
    # for j = index[2]:index[2]+2r
        # @simd for i = index[1]:index[1]+2r 
            # sum += data.modeldata[1].extended[i, j]
        # end
    # end
    # sum
end

@inline update_cell!(hood::DispersalNeighborhood{:outwards}, model, data, state::AbstractFloat,
                   hood_index, dest_index, layers, args...) = begin
    @inbounds propagules = state * model.fraction * hood.kernel[hood_index...]
    @inbounds data.dest[dest_index...] += propagules
    propagules
end

@inline update_cell!(hood::DispersalNeighborhood{:outwards}, model, data, state::Integer,
                   hood_index, dest_index, layers, args...) = begin
    @inbounds spec_rand(source, Float64, args...) * hood.kernel[hood_index...] > model.prob_threshold || return zero(state)

    @inbounds data.dest[dest_index...] = oneunit(state)
    oneunit(state)
end

@inline update_cell!(hood::DispersalNeighborhood{:inwards}, model, data, state, hood_index,
             source_index, args...) = begin
    @inbounds hood.kernel[hood_index...] * data.source[source_index...]
end

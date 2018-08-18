

"""
    DispersalNeighborhood(; f=d -> exponential(d, 1), radius=3, overflow=Skip())
Constructor for neighborhoods, using a dispersal kernel function and a cell radius.

### Keyword Arguments:
- `f::Function`: any function that accepts a Number argument and returns a propbability between 0.0 and 1.0
- `radius::Integer`: a positive integer
- `overflow = Skip()
"""
build_dispersal_kernel(f, param, init, cellsize, r) = begin
    sze = 2r + 1
    kernel = similar(init, Float64, sze, sze)
    # Paper: l. 97
    for y = -r:r, x = -r:r
        kernel[y+r+1, x+r+1] = f(sqrt(y^2 + x^2) * cellsize, param)
    end
    # Normalise
    kernel ./= sum(kernel)
end

# Paper: l. 96
exponential(d, a) = exp(-d / a)

HudginsDispersalGrid(init, suit, human) = begin
    kernel = similar(init, Matrix{Float64})
    h, w = size(init)
    for i in 1:h, j in 1:w t = similar(init, Float64)
        d = similar(init, Float64)
        f = similar(init, Float64)
        for ii in 1:h, jj in 1:w
            d[ii, jj] = sqrt((i - ii)^2 + (j - jj)^2)
            ZI = -0.8438 * suit[ii, jj] + -0.1378 * human[ii, jj]
            f[ii, jj] = 2 * 1.1248 * exp(ZI)/(1+exp(ZI))
        end
        for ii in 1:h, jj in 1:w
            t[ii, jj] = exp(-d[ii, jj] * f[ii,jj])/sum(exp.(-d .* f))
        end
        kernel[i, j] = t
    end
    HudginsDispersalGrid(kernel)
end

"""
Dispersal function taken from Hudgins, 'Predicting the spread of all
invasive forest pests in the United States', 2017
"""
neighborhood(hood::HudginsDispersalGrid, model::HudginsDispersal, state, index, t, source, dest, args...) = begin
    # Ignore cells below the population threshold
    state > output.pop_threshold || return zero(eltype(kernel[1,1]))

    # Setup
    height, width = size(hood.kernel)
    propagules = zero(eltype(kernel[1,1]))

    # Disperse to the entire grid
    for i = 1:height, j = 1:width
        # Skip the state cell
        i, j == index && continue
        # Invade the cell
        dest[i, j] += kernel[index...][i, j]
        propagules += kernel[index...][i, j]
    end
    propagules
end

"""
    neighbors(hood::DispersalNeighborhood, state, index, t, source, dest, args...)
Returns nieghbors for a [`DispersalNeighborhood`](@ref), looping over
the array of dispersal propabilities.
"""
neighbors(hood::DispersalNeighborhood, model, state, index, t, source, dest, args...) = begin
    cc = 0.0
    r = hood.radius

    # Loop over dispersal kernel grid dimensions
    # TODO use an OffsetArray with 0,0 as the center
    for b = -r:r, a = -r:r
        # Ignore the current center cell
        b == 0 && a == 0 && continue
        # Check boundaries
        y, x, is_inbounds = inbounds(index .+ (b, a), size(source), hood.overflow)
        is_inbounds || continue
        # Update cumulative value, and cell value for outwards dispersal
        cc += update_cell!(hood, model, state, t, source, dest, (b, a) .+ (r + 1), (y, x), args...)
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
                   t, source, dest, hood_index, dest_index, args...) = begin
    suitability(model.layers, dest_index, t) > model.suitability_threshold || return zero(state)
    rand() * hood.kernel[hood_index...] > model.prob_threshold || return zero(state)
    # Invade the cell
    dest[dest_index...] = oneunit(state)
end

update_cell!(hood::DispersalNeighborhood{:inwards}, model, state, t, source, dest, hood_index, source_index, args...) = 
    source[source_index...] * hood.kernel[hood_index...] 


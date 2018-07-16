
@mix @with_kw struct Dispersal{L,S,T,N}
    # "[`AbstractLayers`](@ref) or a single [AbstractLayer](@ref). The default is `nothing`."
    layers::L = nothing
    # "A number or Unitful.jl distance."
    cellsize::S = 1.0
    # "Minimum habitat suitability index."
    suitability_threshold::T = 0.1
    # "Neighborhood to disperse to or from"
    neighborhood::N = DispersalNeighborhood(cellsize=cellsize)
end

@mix struct Probabilistic{P}
    # "A real number between one and zero."
    prob_threshold::P = 0.1
end

"Extend to modify [`InwardsLocalDispersal`](@ref)"
abstract type AbstractInwardsDispersal <: AbstractModel end

"""
Local dispersal within a [`DispersalNeighborhood`](@ref) or other neighborhoods.
Inwards dispersal calculates dispersal *to* the current cell from cells in the neighborhood.
"""
@Probabilistic @Dispersal struct InwardsLocalDispersal{} <: AbstractInwardsDispersal end

"Extend to modify [`OutwardsLocalDispersal`](@ref)"
abstract type AbstractOutwardsDispersal <: AbstractPartialModel end

"""
Local dispersal within a [`DispersalNeighborhood`](@ref)

Outwards dispersal calculates dispersal *from* the current cell to cells
in its neighborhood. This should be more efficient than inwards
dispersal when a small number of cells are occupied, but less efficient when a large
proportion of the grid is occupied.
"""
@Probabilistic @Dispersal struct OutwardsLocalDispersal{} <: AbstractOutwardsDispersal end

@Dispersal struct HudginsDispersal{} <: AbstractOutwardsDispersal
    pop_threshold::Float64 = 0.0006227
    growthrate::Float64 = 2.4321
end

"Extend to modify [`JumpDispersal`](@ref)"
abstract type AbstractJumpDispersal <: AbstractPartialModel end

"Jump dispersal within a [`DispersalNeighborhood`](@ref)] or other neighborhoods."
@Probabilistic @Dispersal struct JumpDispersal{L,S} <: AbstractJumpDispersal
    "A number or Unitful.jl distance with the same units as cellsize"
    spotrange::S = 30.0
end

"Inherit to extend human dispersal."
abstract type AbstractHumanDispersal <: AbstractPartialModel end
"Human dispersal model."
@Probabilistic @Dispersal struct HumanDispersal{L,S} <: AbstractHumanDispersal 
    # "A number or Unitful.jl distance."
    spotrange::Int = 30.0
end


"Neighborhoods for dispersal"
abstract type AbstractDispersalNeighborhood <: AbstractNeighborhood end

"""
A neighborhood built from a dispersal kernel function and a cell radius.
Can be built directly by passing in the array, radius and overflow
arguments, but preferably use the keyword constructor to build the array from
a dispersal kernel function.
"""
struct DispersalNeighborhood{K,S} <: AbstractDispersalNeighborhood
    kernel::K
    radius::Int
    overflow::S
end

"""
    DispersalNeighborhood(; f=d -> exponential(d, 1), radius=3, overflow=Skip())
Constructor for neighborhoods, using a dispersal kernel function and a cell radius.

### Keyword Arguments:
- `f::Function`: any function that accepts a Number argument and returns a propbability between 0.0 and 1.0
- `radius::Integer`: a positive integer
- `overflow = Skip()
"""
DispersalNeighborhood(; f=d -> exponential(d, 1), radius=3, cellsize=1, overflow=Skip()) = begin
    r = radius
    size = 2r + 1
    kernel = zeros(Float64, size, size)
    for i = -r:r, j = -r:r
        kernel[i+r+1, j+r+1] = f(norm([i, j]) * cellsize)
    end
    # Zero central cell
    kernel[r + 1, r + 1] = 0.0
    # Normalise
    kernel ./= sum(kernel)
    DispersalNeighborhood(kernel, radius, overflow)
end

exponential(d, a) = e^-d * a

abstract type AbstractDispersalGrid <: AbstractDispersalNeighborhood end

struct HudginsDispersalGrid{K} <: AbstractDispersalGrid
    kernel::K
end

HudginsDispersalGrid(init, suit, human) = begin
    kernel = similar(init, Matrix{Float64})

    h, w = size(init)
    for i in 1:h, j in 1:w
        t = similar(init, Float64)
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

# more dipersal kernel functions here

"""
    pressure(model, cc)
Calculates the propagule pressure from the output of a neighborhood.
"""
pressure(model, cc) = rand()^model.prob_threshold > (1 - cc) / 1

"""
    rule(model::AbstractInwardsDispersal, state, index, t, args...)
Runs rule for of [`AbstractInwardsDispersal`](@ref) dispersal.

The current cell is invaded if there is pressure from surrounding cells and
suitable habitat. Otherwise it keeps its current state.
"""
rule(model::AbstractInwardsDispersal, state, index, t, args...) = begin
    # Exit unless cell habitat is suitabile for invasion
    suit = suitability(model.layers, index, t)
    suit >= model.suitability_threshold || return zero(state)

    # Combine neighborhood cells into a single scalar
    cc = neighbors(model.neighborhood, state, index, t, args...)

    # Set to occupied if suitable habitat and enough pressure from neighbors
    pressure(model, cc) ? 1 : state
end

"""
    rule(model::AbstractOutwardsDispersal, state, index, t, source, dest, args...)
Runs rule for of [`AbstractOutwardsDispersal`](@ref) dispersal.

Surrounding cells are invaded if the current cell is occupied and they have
suitable habitat. Otherwise they keeps their current state.
"""
rule(model::AbstractOutwardsDispersal, state::Integer, index, t, source, dest, args...) = begin
    state == zero(state) && return # Ignore empty cells 

    invade(model.neighborhood, model, state, index, t, source, dest, args...)

    # Set dest cell state to occupied
    dest[index...] = oneunit(state)
end

"""
Dispersal function taken from Hudgins, 'Predicting the spread of all
invasive forest pests in the United States', 2017
"""
rule(model::AbstractOutwardsDispersal, state, index, t, source, dest, args...) = begin
    state == zero(state) && return # Ignore empty cells
    # Grow population - easier to do at the start than the end
    state *= model.growthrate

    propagules = invade(model.neighborhood, model, state, index, t, source, dest, args...)

    # Write the new popuation size to the dest array
    dest[row, col] = state - propagules
end

invade(hood::DispersalNeighborhood, model, state, index, t, source, dest, args...) = begin
    r = hood.radius
    propagules = 0

    # Disperse to the local neighborhood
    for a = -r:r, b = -r:r
        # Skip the state cell
        a == 0 && b == 0 && continue
        # Coordinates are on the grid
        i, j, is_inbounds = inbounds(index .+ (a, b), size(dest), hood.overflow)
        is_inbounds || continue
        # Ignore already populated c
        dest[i, j] == 1 && continue
        # Habitat is suitabile
        suitability(model.layers, (i, j), t) > model.suitability_threshold || continue
        # Invasion pressure is above the threshold
        rand() * hood.kernel[a + r + 1, b + r + 1] > model.prob_threshold || continue 

        # Invade the cell
        dest[i, j] = 1
        propagules += 1
    end
    propagules
end

invade(hood::HudginsDispersalGrid, model::HudginsDispersal, state, index, t, source, dest, args...) = begin
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
    rule(model::AbstractJumpDispersal, state, index, t, source, dest, args...)
Long range rule for [`AbstractJumpDispersal`](@ref). A random cell
within the spotrange is invaded if it is suitable.
"""
rule(model::AbstractJumpDispersal, state, index, t, source, dest, args...) = begin
    # Ignore empty cells
    state > zero(state) || return
    # Random dispersal events
    rand() < model.prob_threshold || return

    # Calculate maximum spotting distance
    range = -model.spotrange:model.spotrange ./ model.cellsize
    # Randomly select actual spotting distance
    spot = tuple(round.(Int, rand(range, 2) .+ index)...)

    # Update spotted cell if it's on the grid and suitable habitat
    row, col, is_inbounds = inbounds(spot, size(dest), Skip())
    if is_inbounds && suitability(model.layers, (row, col), t) > model.suitability_threshold
        dest[row, col] = oneunit(state)
    end
end

"""
    rule(model::AbstractHumanDispersal, state, index, t, source, dest, args...)
Simulates human dispersal, weighting dispersal probability based on human
population in the source cell.
"""
rule(model::AbstractHumanDispersal, state, index, t, source, dest, args...) = begin
    # Ignore empty cells
    state > zero(state) || return

    rand() < model.prob_threshold * human_impact(model.layers, index, t) || return

    # Calculate maximum spotting distance
    range = -model.spotrange:model.spotrange ./ model.cellsize
    # Randomly select actual spotting distance
    spot = tuple(round.(Int, rand(range, 2) .+ index)...)

    # Update spotted cell if it's on the grid and suitable habitat
    row, col, is_inbounds = inbounds(spot, size(dest), Skip())
    if is_inbounds && suitability(model.layers, (row, col), t) * human_impact(model.layers, (row, col), t) > model.suitability_threshold
        dest[row, col] = oneunit(state)
    end
end


"""
    neighbors(hood::DispersalNeighborhood, state, index, t, source, dest, args...)
Returns nieghbors for a [`DispersalNeighborhood`](@ref), looping over
the array of dispersal propabilities.
"""
neighbors(hood::DispersalNeighborhood, state, index, t, source, dest, args...) = begin
    # Setup
    height, width = size(source)
    row, col = index
    cc = 0.0
    radius = hood.radius

    # Loop over dispersal kernel grid dimensions
    for b = -radius:radius, a = -radius:radius
        # Ignore the current center cell
        a == 0 && b == 0 && continue
        # Check boundaries
        p, q, is_inbounds = inbounds((row + b, col + a), (height, width), hood.overflow)
        is_inbounds || continue
        # Update cumulative value with cell value and dispersal kernel
        cc += source[p, q] * hood.kernel[b + radius + 1, a + radius + 1]
    end
    return cc
end

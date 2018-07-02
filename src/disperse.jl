
@mix @with_kw struct MixLocalDispersal{N,L,P,S} 
    neighborhood::N = DispersalNeighborhood()
    "A real number between one and zero."
    prob::P = 0.1
    "[AbstractLayers](@ref) or a single [`AbstractLayer`](@ref)."
    layers::L = nothing
    "A number or Unitful.jl distance."
    cellsize::S = 1.0
end

"""
Extend to modify [`LocalInwardsDispersal`](@ref)
"""
abstract type AbstractInwardsLocalDispersal <: AbstractModel end

"""
Local dispersal within a [`DispersalNeighborhood`](@ref) or other [AbstractNeighborhood](@ref).
Inwards dispersal calculates dispersal *to* the current cell from cells in the neighborhood.
"""
@MixLocalDispersal struct InwardsLocalDispersal{} <: AbstractInwardsLocalDispersal end

"""
Extend to modify [`LocalOutwardsDispersal`](@ref)
"""
abstract type AbstractOutwardsLocalDispersal <: AbstractPartialModel end

"""
Local dispersal within a [`DispersalNeighborhood`](@ref)

Outwards dispersal calculates dispersal *from* the current cell to cells 
in its neighborhood. This should be more efficient than inwards 
dispersal when a small number of cells are occupied, but less efficient when a large 
proportion of the grid is occupied.
"""
@MixLocalDispersal struct OutwardsLocalDispersal{} <: AbstractOutwardsLocalDispersal end

"""
Extend to modify [`JumpDispersal`](@ref)
"""
abstract type AbstractJumpDispersal <: AbstractPartialModel end

"""
Jump dispersal within a [`DispersalNeighborhood`](@ref)] or other [AbstractNeighborhood](@ref)

"""
@with_kw struct JumpDispersal{L,S} <: AbstractJumpDispersal
    "A number or Unitful.jl distance."
    spotrange::Int = 30.0
    "A real number between one and zero."
    prob::Float64 = 0.01
    "[`AbstractLayers`](@ref) or a single [AbstractLayer](@ref). The default is `nothing`."
    layers::L = nothing
    "A number or Unitful.jl distance."
    cellsize::S = 1.0
end

"""
Inherit to exten human dispersal
"""
abstract type AbstractHumanDispersal <: AbstractPartialModel end
"""
Human dispersal model.
"""
@with_kw struct HumanDispersal{L,S} <: AbstractHumanDispersal
    human::Float64
    layers::L
    cellsize::S = 1.0
end


"""
Neighborhoods for dispersal
"""
abstract type AbstractDispersalNeighborhood <: AbstractNeighborhood end

"""
A neighborhood built from a dispersal kernel function and a cell radius.
Can be built directly by passing in the array, radius and overflow
arguments, but preferably use the keyword constructor to build the array from 
a dispersal kernel function.

### Arguments:
- `dispkernel::AbstractArray{T,2}`
- `overflow::AbstractOverflow`
"""
struct DispersalNeighborhood{K,S} <: AbstractDispersalNeighborhood
    dispkernel::K
    radius::Int
    overflow::S
end

"""
    DispersalNeighborhood(; f=d -> exponential(d, 1), radius=3, overflow=Skip())
Constructor for neighborhoods, using a dispersal kernel function and a cell radius.

### Keyword Arguments:
- `f::Function`: any function that accepts a Number argument and returns a propbability between 0.0 and 1.0
- `radius::Integer`: a positive integer
"""
DispersalNeighborhood(; f=d -> exponential(d, 1), radius=3, overflow=Skip()) = begin
    dispkernel = build_dispersal_kernel(f, radius)
    DispersalNeighborhood(dispkernel, radius, overflow)
end

"""
    build_dispersal_kernel(f, r)
Accepts a dispersal kernel function and integer radius,
and returns an array of probabilities, of size  2r + 1 * 2r + 1.
"""
build_dispersal_kernel(f, r) = begin
    size = 2r + 1
    grid = zeros(Float64, size, size)
    for i = -r:r, j = -r:r
        d = sqrt(i^2 + j^2)
        grid[i + r + 1, j + r + 1] = f(d)
    end
    grid
end

exponential(d, a) = e^-d * a

# more dipersal kernel functions here

"""
    pressure(model, cc)
Calculates the propagule pressure from the output of a neighborhood.
"""
pressure(model, cc) = rand()^model.prob > (1 - cc) / 1

"""
    rule(model::AbstractInwardsLocalDispersal, state, index, t, args...)
Runs rule for of [`AbstractLocalInwardsDispersal`](@ref) dispersal. 

The current cell is invaded if there is pressure from surrounding cells and 
suitable habitat. Otherwise it keeps its current state.
"""
rule(model::AbstractInwardsLocalDispersal, state, index, t, args...) = begin
    # Exit unless cell habitat is suitabile for invasion
    suitability(model.layers, index..., t) > 0.1 || return state

    # Combine neighborhood cells into a single scalar
    cc = neighbors(model.neighborhood, state, index, t, args...)

    # Set to occupied if suitable habitat and enough pressure from neighbors
    pressure(model, cc) ? 1 : state
end

"""
    rule(model::AbstractOutwardsLocalDispersal, state, index, t, source, dest, args...)
Runs rule for of [`AbstractLocalOutwardsDispersal`](@ref) dispersal. 

Surrounding cells are invaded if the current cell is occupied and they have 
suitable habitat. Otherwise they keeps their current state.
"""
rule(model::AbstractOutwardsLocalDispersal, state, index, t, source, dest, args...) = begin
    # Ignore empty cells completely
    state == zero(state) && return

    # Set dest cell state to occupied
    dest[index...] = oneunit(state)

    # Setup
    height, width = size(dest)
    row, col = index
    hood = model.neighborhood

    # Disperse to surrounding cells
    for a = -hood.radius:hood.radius, b = -hood.radius:hood.radius
        a == 0 && b == 0 && continue
        p, q, is_inbounds = inbounds((row + b, col + a), (height, width), hood.overflow)
        # TODO 0.1 as a parameter. Is it best described as the 'threshold'?
        if is_inbounds && rand() > model.prob && suitability(model.layers, p, q, t) > 0.1
            dest[p, q] |= 1
        end
    end
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
    rand() < model.prob || return

    # Calculate maximum spotting distance
    range = -model.spotrange:model.spotrange ./ model.cellsize
    # Randomly select actual spotting distance
    spot = tuple(round.(Int, rand(range, 2) .+ index)...)

    # Update spotted cell if it's on the grid and suitable habitat
    row, col, is_inbounds = inbounds(spot, size(dest), Skip())
    if is_inbounds && suitability(model.layers, row, col, t) > 0.1 # TODO 0.1 as a parameter - 'threshold'?
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

    rand() < model.prob * human_impact(model.layers, row, col, t) || return

    # Calculate maximum spotting distance
    range = -model.spotrange:model.spotrange ./ model.cellsize
    # Randomly select actual spotting distance
    spot = tuple(round.(Int, rand(range, 2) .+ index)...)

    # Update spotted cell if it's on the grid and suitable habitat
    row, col, is_inbounds = inbounds(spot, size(dest), Skip())
    if is_inbounds && suitability(model.layers, row, col, t) > 0.1
        source[row, col] = oneunit(state)
    end
end


"""
    neighbors(hood::DispersalNeighborhood, state, index, t, source, dest, args...) = begin
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
    for a = -radius:radius, b = -radius:radius
        # Ignore the current center cell
        a == 0 && b == 0 && continue
        # Check boundaries
        p, q, is_inbounds = inbounds((row + b, col + a), (height, width), hood.overflow)
        is_inbounds || continue
        # Update cumulative value with cell value and dispersal kernel
        cc += source[p, q] * hood.dispkernel[a + radius + 1, b + radius + 1]
    end
    return cc
end

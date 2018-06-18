"""
$(TYPEDEF)
Extend to modify [`LocalDispersal`](@ref)
"""
abstract type AbstractLocalDispersal <: AbstractModel end

"""
$(TYPEDEF)

Local dispersal within a [`DispersalNeighborhood`](@ref)] or other [AbstractNeighborhood](@ref)

### Keyword Arguments
- neighborhood: = DispersalNeighborhood()
- prob: A real number between one and zero.
- layers: [AbstractLayers](@ref) or a single [`AbstractLayer`](@ref). The default is `nothing`.
- cellsize::S = A number or Unitful.jl distance.
"""
@with_kw struct LocalDispersal{N,L,P,S} <: AbstractLocalDispersal
    neighborhood::N = DispersalNeighborhood()
    prob::P = 0.1
    layers::L = nothing
    cellsize::S = 1.0
end

"""
$(TYPEDEF)
Extend to modify [`JumpDispersal`](@ref)
"""

abstract type AbstractJumpDispersal <: AbstractInPlaceModel end
"""
$(TYPEDEF)
Local dispersal within a [`DispersalNeighborhood`](@ref)] or other [AbstractNeighborhood](@ref)

### Keyword Arguments
- spotrange: A number or Unitful.jl distance.
- prob: A real number between one and zero.
- layers: [`AbstractLayers`](@ref) or a single [AbstractLayer](@ref). The default is `nothing`.
- cellsize::S = A number or Unitful.jl distance.
"""
@with_kw struct JumpDispersal{L,S} <: AbstractJumpDispersal
    spotrange::Int = 30.0
    prob::Float64 = 0.01
    layers::L = nothing
    cellsize::S = 1.0
end

"""
$(TYPEDEF)
Inherit to exten human dispersal
"""
abstract type AbstractHumanDispersal <: AbstractInPlaceModel end
"""
$(TYPEDEF)
Human dispersal model.
$(FIELDS)
"""
@with_kw struct HumanDispersal{L,S} <: AbstractHumanDispersal
    human::Float64
    layers::L
    cellsize::S = 1.0
end


"""
$(TYPEDEF)
Neighborhoods for dispersal
"""
abstract type AbstractDispersalNeighborhood <: AbstractNeighborhood end

"""
$(TYPEDEF)
A neighborhood built from a dispersal kernel function and a cell radius.

### Arguments:
- `dispkernel::AbstractArray{T,2}`
- `overflow::AbstractOverflow`
"""
struct DispersalNeighborhood{K,S} <: AbstractDispersalNeighborhood
    dispkernel::K
    overflow::S
end

"""
Build a neighborhood from a dispersal kernel function and a cell radius.

### Keyword Arguments:
- `f::Function`: any function that accepts a Number argument and returns a propbability between 0.0 and 1.0
- `radius::Integer`: a positive integer 
"""
DispersalNeighborhood(; f = d -> exponential(d, 1), radius = 3, 
                      overflow = Skip()) = begin
    dispkernel = build_dispersal_kernel(f, radius) 
    DispersalNeighborhood(dispkernel, overflow)
end

"""
$(SIGNATURES)
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
$(METHODLIST)
Calculates the propagule pressure from the output of a neighborhood.
"""
function pressure end

pressure(model::AbstractLocalDispersal, cc) = rand()^model.prob > (1 - cc) / 1

"""
$(SIGNATURES)
Short range rule for [`AbstractLocalDispersal`](@ref) dispersal. Cells are invaded 
if there is pressure and suitable habitat, otherwise left as-is.
"""
rule(model::AbstractLocalDispersal, state, index, t, args...) = begin
    cc = neighbors(model.neighborhood, state, index, t, args...)
    suitable = suitability(model.layers, index..., t) > 0.1
    suitable && pressure(model, cc) ? oneunit(state) : state
end

"""
$(SIGNATURES)
Long range rule for [`AbstractJumpDispersal`](@ref). A random cell
within the spotrange is invaded if it is suitable.
"""
rule(model::AbstractJumpDispersal, state, index, t, source, args...) = begin
    if state > zero(state) && rand() < model.prob
        range = -model.spotrange:model.spotrange ./ model.cellsize
        spot = tuple(round.(Int, rand(range, 2) .+ index)...)
        row, col, ok = inbounds(spot, size(source), Skip())
        if ok && suitability(model.layers, row, col, t) > 0.1
            source[row, col] = oneunit(state)
        end
    end
end

rule(f::AbstractHumanDispersal, model, state, index, t, source, args...) = begin
    if state > zero(state)
        if ok && suitability(model.layers, row, col, t) 
            source[row, col] = oneunit(state)
        end
    end
end

"""
$(SIGNATURES)
Returns nieghbors for a [`DispersalNeighborhood`](@ref), looping over
the array of dispersal propabilities.
"""
neighbors(h::DispersalNeighborhood, state, index, t, source, args...) = begin
    height, width = size(source)
    row, col = index
    r = div(size(h.dispkernel, 1) - 1, 2)
    cc = 0.0
    # loop over dispersal kernel grid dimensions
    for a = -r:r 
        for b = -r:r
            # ignore the current center cell
            a == 0 && b == 0 && continue
            p, q, inb = inbounds((row + b, col + a), (height, width), h.overflow)
            inb || continue
            cc += source[p, q] * h.dispkernel[a + r + 1, b + r + 1]
        end
    end
    return cc
end

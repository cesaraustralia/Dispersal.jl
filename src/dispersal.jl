
abstract type AbstractLocalDispersal <: AbstractCellular end
@with_kw struct LocalDispersal{N,L,S} <: AbstractLocalDispersal
    neighborhood::N = DispersalNeighborhood()
    prob::Float64 = 0.9
    layers::L
    cellsize::S = 1.0
end

abstract type AbstractJumpDispersal <: AbstractInPlaceCellular end
@with_kw struct JumpDispersal{L,S} <: AbstractJumpDispersal
    prob::Float64 = 0.01
    spotrange::Int = 30.0
    layers::L
    cellsize::S = 1.0
end

abstract type AbstractHumanDispersal <: AbstractInPlaceCellular end
@with_kw struct HumanDispersal{L,S} <: AbstractHumanDispersal
    human::Float64
    layers::L
    cellsize::S = 1.0
end


abstract type AbstractDispersalNeighborhood <: AbstractNeighborhood end
struct DispersalNeighborhood{K,S} <: AbstractDispersalNeighborhood
    dispkernel::K
    overflow::S
end

DispersalNeighborhood(; f = d -> exponential(d, 1), radius = 3, 
                      overflow = Skip()) = begin
    dispkernel = build_dispersal_kernel(f, radius) 
    DispersalNeighborhood(dispkernel, overflow)
end

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


pressure(::AbstractLocalDispersal, cc) = rand() > (1 - cc) / 1

"""
Short range rule for dispersal kernels. Cells are invaded if there is pressure and 
suitable habitat, otherwise left as-is.
"""
rule(model::LocalDispersal, state, index, t, args...) = begin
    cc = neighbors(model.neighborhood, state, index, t, args...)
    press = pressure(model, cc) 
    suitable = suitability(model.layers, index..., t) > 0.1
    press && suitable ? oneunit(state) : state
end

"""
Long range rule for dispersal kernels. Cells are invaded if there is pressure and 
suitable habitat, otherwise left as-is.
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


"Extend to modify [`HumanDispersal`](@ref)"
abstract type AbstractHumanDispersal <: AbstractPartialModel end

# "Human dispersal model."
@Probabilistic struct HumanDispersal{PC} <: AbstractHumanDispersal 
    precalc::PC = [] | false | _
end

HumanDispersal(precalc::PC, cellsize::CS, prob_threshold::PT) where {PC,CS,PT} =
    HumanDispersal{PC,CS,PT}(precalc, cellsize, prob_threshold)

abstract type AbstractCell end

mutable struct CellMagnitude{M,I} <: AbstractCell
    magnitude::M
    ind::I
end

mutable struct CellInterval{P,M,I} <: AbstractCell
    cumprop::P
    fraction::P
    magnitude::M
    ind::I
end

import Base: isless, +, zero
isless(x::CellMagnitude, y::CellMagnitude) = isless(x.magnitude, y.magnitude)
isless(x::CellMagnitude, y) = isless(x.magnitude, y)
isless(x, y::CellMagnitude) = isless(x, y.magnitude)
isless(x::CellInterval, y::CellInterval) = isless(x.cumprop, y.cumprop)
isless(x::CellInterval, y) = isless(x.cumprop, y)
isless(x, y::CellInterval) = isless(x, y.cumprop)

zero(x::Type{CellMagnitude{T,Tuple{A,B}}}) where {T,A,B} = zero(x.magnitude)

+(x::CellMagnitude, y::CellMagnitude) = +(x.magnitude, y.magnitude)
+(x, y::CellMagnitude) = +(x, y.magnitude)
+(x::CellMagnitude, y) = +(x.magnitude, y)

" Precalculate a dispersal shortlist for each cell "
precalc_human_dispersal(human::A, cellsize, take) where A <: AbstractArray{T} where T = begin
    h, w = size(human)
    indices = broadcastable_indices(Int32, human)
    dist = distances(human) .* cellsize
    s = similar(human)
    take = min(take, length(human))
    magnitudes = Matrix{CellMagnitude{Float32,Tuple{Int32,Int32}}}(undef, size(human)...)
    broadcast!(index -> CellMagnitude(0.0f0, index), magnitudes, indices)
    flat_magnitudes = Vector{CellMagnitude{Float32,Tuple{Int32,Int32}}}(undef, size(human, 1) * size(human, 2))
    top_magnitudes = Vector{CellMagnitude{Float32,Tuple{Int32,Int32}}}(undef, take)
    intervals = [CellInterval(0.0f0, 0.0f0, 0.0f0, Int32.((0, 0))) for i in 1:take]
    precalc = [Vector{CellInterval{Float32,Float32,Tuple{Int32,Int32}}}(undef, take) for i in 1:size(human, 1), j in 1:size(human, 2)]
    props = similar(human)

    for i in 1:size(human, 1), j in 1:size(human, 2)
        broadcast(build_cell_pop_index, magnitudes, i, j, indices, (human,), (dist,)) # 4
        for n in 1:size(magnitudes, 1) * size(magnitudes, 2)
            flat_magnitudes[n] = magnitudes[n]
        end
        partialsort!(flat_magnitudes, take, rev=true)
        top_magnitudes .= flat_magnitudes[1:take]
        mag_sum::Float32 = sum(top_magnitudes)
        realmag_sum::Float32 = sum(magnitudes)
        cumprop = 0.0f0
        for (n, m) in enumerate(reverse(top_magnitudes))
            prop = m.magnitude / mag_sum
            cumprop += prop
            precalc[i, j][n] = CellInterval(cumprop, prop, m.magnitude, m.ind)
        end
        props[i, j] = mag_sum / realmag_sum
    end

    precalc, props
end

build_cell_pop_index(m, i, j, (ii, jj), human, dist) = begin
    m.magnitude = exponential(dist[abs(i - ii) + 1, abs(j - jj) + 1], 50) * human[ii, jj]
end

""" 
Populate a matrix with the values of a list of cells. 
This lets you view the contents of a cell in an output display.

## Arguments:
`a`: A matrix of the same size the precalculation was performed on
`cells`: A vector of [`CellInterval`](@ref)
"""
populate!(a::AbstractMatrix, cells::AbstractVector{<:CellInterval}) = begin
    for cell in cells 
        a[cell.ind...] = cell.fraction
    end
    a
end
populate!(a::AbstractMatrix{<:Integer}, cells::AbstractVector{<:CellInterval}) = begin
    for cell in cells 
        a[cell.ind...] = 1
    end
    a
end

populate(cells::AbstractVector{<:CellInterval}, sze) = populate!(zeros(sze...), cells)


"""
    rule(model::AbstractHumanDispersal, state, index, t, source, dest, args...)
Simulates human dispersal, weighting dispersal probability based on human
population in the source cell.
"""
rule!(model::AbstractHumanDispersal, data, state, index, args...) = begin
    # Ignore empty cells
    state > zero(state) || return
    # Randomise dispersal
    rand() < model.prob_threshold || return
    
    # Randomly choose a cell to disperse to from the precalculated human dispersal distribution
    shortlist = model.precalc[index...]
    dest_id = min(length(shortlist), searchsortedfirst(shortlist, rand()))
    dest_index = shortlist[dest_id].ind
    # Disperse to the cell
    update_cell!(model, data, state, dest_index)
    state
end

update_cell!(model::AbstractHumanDispersal, data, state::AbstractFloat, dest_index) = 
    data.dest[dest_index...] += oneunit(state)
update_cell!(model::AbstractHumanDispersal, data, state::Integer, dest_index) = 
    data.dest[dest_index...] = oneunit(state)

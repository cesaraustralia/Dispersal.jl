"Extend to modify [`HumanDispersal`](@ref)"
abstract type AbstractHumanDispersal <: AbstractPartialModel end

# "Human dispersal model."
@Probabilistic @SpotRange @Dispersal @Suitability struct HumanDispersal{P} <: AbstractHumanDispersal 
    precalc::P = [] | false | _
end

HumanDispersal(prob_threshold::PT, spotrange::SR, cellsize::CS, precalc::PC) where {PT,SR,CS,PC} = begin
    precalc = precalc_human_dispersal(cellsize, spotrange)
    HumanDispersal{PT,SR,CS,typeof(precalc)}(cellsize, prob_threshold, spotrange, precalc)
end


mutable struct CellMagnitude{M,I}
    magnitude::M
    ind::I
end

mutable struct CellInterval{P,M,I}
    cumprop::P
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

+(x::CellInterval, y::CellInterval) = +(x.prop, y.prop)
+(x, y::CellInterval) = +(x, y.prop)
+(x::CellInterval, y) = +(x.prop, y)
+(x::CellMagnitude, y::CellMagnitude) = +(x.magnitude, y.magnitude)
+(x, y::CellMagnitude) = +(x, y.magnitude)
+(x::CellMagnitude, y) = +(x.magnitude, y)

build_cell_pop_index(m, i, j, ii, jj, human, dist) = begin
    m.magnitude = exponential(dist[abs(i - ii) + 1, abs(j - jj) + 1], 50) * human[ii, jj]
end

precalc_human_dispersal(human::A, cellsize, take) where A <: AbstractArray{T} where T = begin
    h, w = size(human)
    rows, cols = broadcastable_indices(Int32, human)
    dist = distances(human) .* cellsize
    s = similar(human)
    take = min(take, length(human))
    precalc = [Array{CellInterval{Float32,Float32,Tuple{Int32,Int32}}}(take) for i in 1:size(human, 1), j in 1:size(human, 2)]
    magnitudes = Array{CellMagnitude{Float32,Tuple{Int32,Int32}}}(size(human)...)
    broadcast!((r, c) -> CellMagnitude(0.0f0, (r,c)), magnitudes, rows, cols)
    flat_magnitudes = Array{CellMagnitude{Float32,Tuple{Int32,Int32}}}(size(human, 1) * size(human, 2))
    top_magnitudes = Array{CellMagnitude{Float32,Tuple{Int32,Int32}}}(take)
    intervals = [CellInterval(0.0f0, 0.0f0, Int32.((0, 0))) for i in 1:take]
    props = similar(human)

    for i in 1:size(human, 1), j in 1:size(human, 2)
        broadcast(build_cell_pop_index, magnitudes, i, j, rows, cols, (human,), (dist,)) # 4
        for n in 1:size(magnitudes, 1) * size(magnitudes, 2)
            flat_magnitudes[n] = magnitudes[n]
        end
        partialsort!(flat_magnitudes, take, rev=true)
        top_magnitudes .= flat_magnitudes[1:take]
        mag_sum::Float32 = sum(top_magnitudes)
        realmag_sum::Float32 = sum(magnitudes)
        cumprop = 0.0f0
        for (n, m) in enumerate(top_magnitudes)
            cumprop += m.magnitude / mag_sum
            precalc[i, j][n] = CellInterval(cumprop, m.magnitude, m.ind)
        end
        props[i, j] = mag_sum / realmag_sum
    end

    precalc, props
end

populate!(a::AbstractMatrix, cis::AbstractVector{CellInterval{F,F,Tuple{I,I}}}) where {F,I} = 
    for c in cis 
        a[c.ind...] = c.magnitude
    end

populate!(a::AbstractMatrix{Int}, cis::AbstractVector{CellInterval}) = 
    for c in cis 
        a[c.ind...] = 1
    end


"""
    rule(model::AbstractHumanDispersal, state, row, col, t, source, dest, args...)
Simulates human dispersal, weighting dispersal probability based on human
population in the source cell.
"""
rule!(model::AbstractHumanDispersal, state, row, col, t, source, dest, layers, precalc, args...) = begin
    # Ignore empty cells
    state > zero(state) || return
    
    # Randomly choose a cell to disperse to from the precalculated human dispersal distribution
    shortlist = precalc[row, col]
    dest_id = min(length(shortlist), searchsortedfirst(shortlist, spec_rand(source, Float64, args...)))
    # Disperse to the cell
    dest[shortlist[dest_id].ind...] = oneunit(state)
    state
end

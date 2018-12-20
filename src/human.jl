using SharedArrays

"Extend to modify [`HumanDispersal`](@ref)"
abstract type AbstractHumanDispersal <: AbstractPartialModel end

# "Human dispersal model."
@Probabilistic struct HumanDispersal{PC} <: AbstractHumanDispersal
    precalc::PC = [] | false | _
end

HumanDispersal(precalc::PC, cellsize::CS, prob_threshold::PT) where {PC,CS,PT} =
    HumanDispersal{PC,CS,PT}(precalc, cellsize, prob_threshold)

abstract type AbstractCell end

mutable struct CellGravity{M,I} <: AbstractCell
    gravity::M
    index::I
end

struct CellInterval{P,M,I} <: AbstractCell
    cumprop::P
    fraction::P
    gravity::M
    index::I
end

########################################################################################
# Sorting
#
# isless is used to teratively sort lists and search in funcitons like
# searchsortedfirst() and partialsort!
# we define isless on CellGravity.gravity in order to sort margnitudes in order.
# on CellInterval it is defined on cumprop in order to randomly choose from a list
# maintaining proportion when using searchsortedfirst.

import Base: isless, +
isless(x::CellGravity, y::CellGravity) = isless(x.gravity, y.gravity)
isless(x::CellGravity, y) = isless(x.gravity, y)
isless(x, y::CellGravity) = isless(x, y.gravity)

isless(x::CellInterval, y::CellInterval) = isless(x.cumprop, y.cumprop)
isless(x::CellInterval, y) = isless(x.cumprop, y)
isless(x, y::CellInterval) = isless(x, y.cumprop)

# Adding methods for + allows us to use sum() on arrays of CellGravity
+(x::CellGravity, y::CellGravity) = +(x.gravity, y.gravity)
+(x, y::CellGravity) = +(x, y.gravity)
+(x::CellGravity, y) = +(x.gravity, y)


"""
Precalculate a dispersal shortlist for each cell
"""
precalc_human_dispersal(human_pop::AbstractMatrix, cellsize, shortlist_len, human_exponent, dist_exponent) = begin
    # Precalculate exponentiation of human population matrix
    human = human_pop .^ human_exponent
    # Precalculated distances matrix
    dist = (distances(human) .* cellsize) .^ dist_exponent
    dist[1] =  # mean distance from cell centre 
    # Get matrix dimensions
    h, w = size(human)
    # Get indices to broadcast over
    indices = broadcastable_indices(Int32, human)
    s = similar(human)

    # Limit shortlist cells to the total available
    shortlist_len = min(shortlist_len, length(human))

    # Preallocate memory
    gravitys = Matrix{CellGravity{Float32,Tuple{Int32,Int32}}}(undef, size(human)...)
    broadcast!(index -> CellGravity(0.0f0, index), gravitys, indices)
    gravity_vector = Vector{CellGravity{Float32,Tuple{Int32,Int32}}}(undef, size(human, 1) * size(human, 2))
    gravity_shortlist = Vector{CellGravity{Float32,Tuple{Int32,Int32}}}(undef, shortlist_len)
    interval_shortlist = Vector{CellInterval{Float32,Float32,Tuple{Int32,Int32}}}(undef, shortlist_len)
    # Final precalc array to be returned
    precalc = [Vector{CellInterval{Float32,Float32,Tuple{Int32,Int32}}}(undef, shortlist_len) for i in 1:size(human, 1), j in 1:size(human, 2)]
    # Matrix of Proportions
    props = similar(human)

    # Precalculate human dispersal shortlist for every cell in the grid
    for j = 1:size(human, 2)
        println("Precalculating column : ", j)
        for i = 1:size(human, 1)
            # Calculate the gravityfor all cells in the grid
            broadcast(build_gravity_index, gravitys, i, j, indices, (human,), (dist,)) # 4

            # Arrange gravitys in a vector for 1 dimensional sorting
            for n = 1:size(gravitys, 1) * size(gravitys, 2)
                gravity_vector[n] = gravitys[n]
            end
            # Sort the top shortlist_len gravitys in-place, highest first.
            partialsort!(gravity_vector, shortlist_len, rev=true)
            # Copy sorted gravitys to the shortlist
            gravity_shortlist .= gravity_vector[1:shortlist_len]

            # Sum gravitys in the shortlist
            shortlist_sum::Float32 = sum(gravity_shortlist)
            # Sum all gravitys
            total_sum::Float32 = sum(gravitys)

            # Create a list of intervals from the sorted list of gravitys.
            # This will be used to choose to randomly select cells from the 
            # distribution of gravitys
            cumprop = 0.0f0
            for (n, m) = enumerate(reverse(gravity_shortlist))
                # Calculaate proportion of current gravity in the complete shortlist
                prop = m.gravity / shortlist_sum
                # Track cumulative proportion for use with `searchsortedfirst()`
                cumprop += prop
                interval_shortlist[n] = CellInterval(cumprop, prop, m.gravity, m.index)
            end

            # Update output matrix
            precalc[i, j] .= interval_shortlist
            # Update shortlist proportion matrix to check coverage of the distribution
            props[i, j] = shortlist_sum / total_sum
        end
    end
    precalc, props
end

"""
Calculate the gravity of an individual cell relative to the current cell.

This is a combination of the distance and population of the cell.
"""
@inline build_gravity_index(m, i, j, (ii, jj), human, dist) = begin
    m.gravity = (human[i, j] * human[ii, jj]) / (dist[abs(i - ii) + 1, abs(j - jj) + 1])
end

"""
Populate a matrix from a shortlist of cells from one cell in the precalculated matrix
This lets you view the contents of a cell in an AbstractOutput display.

## Arguments:
`a`: A matrix of the same size the precalculation was performed on
`cells`: A vector of [`CellInterval`](@ref)
"""
populate!(a::AbstractMatrix, cells::AbstractVector{<:CellInterval}) = begin
    for cell in cells
        a[cell.index...] = cell.fraction
    end
    a
end
populate!(a::AbstractMatrix{<:Integer}, cells::AbstractVector{<:CellInterval}) = begin
    for cell in cells
        a[cell.index...] = 1
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
    # Randomly disperse only some proportion of the time
    rand() < model.prob_threshold || return

    shortlist = model.precalc[index...]
    # Randomly choose a cell to disperse to from the precalculated human dispersal distribution
    dest_id = searchsortedfirst(shortlist, rand())
    # Get the index of the selelected cell
    dest_index = shortlist[dest_id].index
    # Disperse to the cell
    update_cell!(model, data, state, dest_index)

    state
end

update_cell!(model::AbstractHumanDispersal, data, state, dest_index) =
    data.dest[dest_index...] += oneunit(state)
update_cell!(model::AbstractHumanDispersal, data, state::Bool, dest_index) =
    data.dest[dest_index...] = oneunit(state)

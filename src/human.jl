using Distributed

"Extend to modify [`HumanDispersal`](@ref)"
abstract type AbstractHumanDispersal <: AbstractPartialModel end

# "Human dispersal model."
@limits @flattenable struct HumanDispersal{PC,H,A,D} <: AbstractHumanDispersal
    precalc::PC              | false | _
    human::H                 | false | _
    par_a::A                 | true  | (0.0, 1.0)
    human_dispersal_probs::D | false | _
end

HumanDispersal(precalc, human, par_a = 1e-5) = begin
    human_dispersal_probs = precalc_human_dispersal_probs(human, par_a)
    HumanDispersal(precalc, human, par_a, human_dispersal_probs)
end

precalc_human_dispersal_probs(human, par_a) = 1 .- 1 ./ exp.(par_a .* human)

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

const Interval = CellInterval{Float32,Float32,Tuple{Int32,Int32}}
const Gravity = CellGravity{Float32,Tuple{Int32,Int32}}
const Index = Tuple{Int64,Int64}

const jobs = RemoteChannel(()->Channel(10))
const results = RemoteChannel(()->Channel{Tuple{Index,Vector{Vector{Interval}},Vector{Float32}}}(10))

"""
Precalculate a dispersal shortlist for each cell
"""
precalc_human_dispersal(human_pop::AbstractMatrix, cellsize, shortlist_len, human_exponent, dist_exponent) = begin
    # Precalculate exponentiation of human population matrix
    human = human_pop .^ human_exponent
    # Precalculated distances matrix
    dist = (distances(human) .* cellsize) .^ dist_exponent
    dist[1] = cellsize/6 * (sqrt(2) + log(1 + sqrt(2))) # mean distance from cell centre
    # Get matrix dimensions
    h, w = size(human)
    # Get indices to broadcast over
    indices = broadcastable_indices(Int32, human)

    # Limit shortlist cells to the total available
    shortlist_len = min(shortlist_len, length(human))

    # Preallocate memory
    gravities = Matrix{Gravity}(undef, size(human)...)
    broadcast!(index -> Gravity(0.0f0, index), gravities, indices)
    gravity_vector = Vector{Gravity}(undef, h * w)
    gravity_shortlist = Vector{Gravity}(undef, shortlist_len)
    interval_shortlist = Vector{Interval}(undef, shortlist_len)
    precalc_col = [Vector{Interval}(undef, shortlist_len) for i in 1:h]
    prop_col = [0.0f0 for i in 1:h]
    data = shortlist_len, indices, human, dist, gravities, gravity_vector, gravity_shortlist, interval_shortlist, precalc_col, prop_col

    # Final precalc array to be returned
    precalcs = [Vector{Interval}(undef, shortlist_len) for i in 1:h, j in 1:w]
    # Matrix of Proportions
    props = similar(human)


    @async for j = 1:w
        put!(jobs, (j, data))
    end

    for p in workers()
        remote_do(do_work, p, jobs, results)
    end

    n = 1
    while n <= w
        j, precalc_col, prop_col = take!(results)
        for i = 1:h
            precalcs[i, j] .= precalc_col[i]
            props[i, j] = prop_col[i]
        end
        n += 1
    end

    precalcs, props
end


function do_work(jobs, results) # define work function everywhere
    while true
        j, data = take!(jobs)
        put!(results, build_cell_precalc(j, data))
    end
end


# Precalculate human dispersal shortlist for every cell in the grid
function build_cell_precalc(j, data)
    shortlist_len, indices, human, dist, gravities, gravity_vector, gravity_shortlist, interval_shortlist, precalc_col, prop_col = data
    for i = 1:size(human, 1)
        # Calculate the gravityfor all cells in the grid
        broadcast(build_gravity_index, gravities, i, j, indices, (human,), (dist,)) # 4

        # Arrange gravities in a vector for 1 dimensional sorting
        for n = 1:size(gravities, 1) * size(gravities, 2)
            gravity_vector[n] = gravities[n]
        end
        # Sort the top shortlist_len gravities in-place, highest first.
        partialsort!(gravity_vector, shortlist_len, rev=true)
        # Copy sorted gravities to the shortlist
        gravity_shortlist .= gravity_vector[1:shortlist_len]

        # Sum gravities in the shortlist
        shortlist_sum::Float32 = sum(gravity_shortlist)
        # Sum all gravities
        total_sum::Float32 = sum(gravities)

        # Create a list of intervals from the sorted list of gravities.
        # This will be used to choose to randomly select cells from the 
        # distribution of gravities
        cumprop = 0.0f0
        for (n, m) = enumerate(reverse(gravity_shortlist))
            # Calculaate proportion of current gravity in the complete shortlist
            prop = m.gravity / shortlist_sum
            # Track cumulative proportion for use with `searchsortedfirst()`
            cumprop += prop
            interval_shortlist[n] = CellInterval(cumprop, prop, m.gravity, m.index)
        end
        prop = shortlist_sum / total_sum
        # Update output matrix
        precalc_col[i] .= interval_shortlist
        # Update shortlist proportion matrix to check coverage of the distribution
        prop_col[i] = prop
    end

    j, precalc_col, prop_col
end

# """
# Calculate the gravity of an individual cell relative to the current cell.

# This is a combination of the distance and population of the cell.
# """

build_gravity_index(m, i, j, (ii, jj), human, dist) = begin
    @inbounds m.gravity = (human[i, j] * human[ii, jj]) / (dist[abs(i - ii) + 1, abs(j - jj) + 1])
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

    dispersalprob = model.human_dispersal_probs[index...]
    meandispersers = round(dispersalprob * state)

    dispersers = meandispersers # deterministic
    # dispersers = Rand.Poisson(meandispersers) # random

    for i = 1:dispersers
        # Randomly choose a cell to disperse to from the precalculated human dispersal distribution
        shortlist = model.precalc[index...]
        dest_id = min(length(shortlist), searchsortedfirst(shortlist, rand()))
        dest_index = shortlist[dest_id].index
        # Disperse to the cell
        update_cell!(model, data, state, dest_index)
    end
    # TODO make method for boolean and float
    data.dest[index...] -= dispersers
    data.dest[index...]
end

update_cell!(model::AbstractHumanDispersal, data, state, dest_index) =
    data.dest[dest_index...] += oneunit(state)
update_cell!(model::AbstractHumanDispersal, data, state::Bool, dest_index) =
    data.dest[dest_index...] = oneunit(state)

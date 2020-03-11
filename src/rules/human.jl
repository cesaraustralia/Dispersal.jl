
"""
CellGravity allows sorting on gravity while keeping records of the original cell coordinate
"""
mutable struct CellGravity{M,I}
    gravity::M
    index::I
end

"""
CellGravity allows ordering a list by the cumulative proportion of the total gravity,
and plotting based on the fraction of total gravity.
"""
struct CellInterval{P,I}
    cumprop::P
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


# Define types used in the precalculation
const Index = Tuple{Int16,Int16}
const Interval = CellInterval{Float32,Index}
const Gravity = CellGravity{Float32,Index}
const Precalc = Union{Vector{Interval},Missing}
const Prop = Union{Float32,Missing}

"""
HumanDispersal Rules human-driven dispersal patterns using population data.

Transport connections between grid cells are calculated using distance and human population,
modified with the `human_exponent` and `dist_exponent` parameters. A shortlist of the most
connected cells is selected for use in the simulation.

The time taken for precalulation will depend on the `scale` argument. Values above 1
will downsample the grid to improve precalulation time and runtime performance. A high
scale value is good for use in a live interface.

## Arguments

$(FIELDDOCTABLE)
"""
@description @limits @flattenable struct HumanDispersal{R,W,HP,CS,S,AG,HE,DE,EA,MD,SL,TS,PC,PR,DP,B} <: PartialRule{R,W}
    # Field                | Flatten | Limits
    human_pop::HP          | false   | _               | _
    cellsize::CS           | false   | _               | _
    scale::S               | false   | _               | _
    aggregator::AG         | false   | _               | "A function that aggregates scaled down cells"
    human_exponent::HE     | true    | (1.0, 3.0)      | "Human population exponent"
    dist_exponent::DE      | true    | (1.0, 3.0)      | "Distance exponent"
    dispersalperpop::EA    | true    | (0.0, 1e-8)     | "Scales the number of dispersing individuals by human activity (ie population^human_exponent)"
    max_dispersers::MD     | true    | (50.0, 10000.0) | "Maximum number of dispersers in a dispersal events"
    shortlist_len::SL      | false   | _               | "Length of dest shortlist"
    timestep::TS           | false   | _               | _
    dest_shortlists::PC    | false   | _               | _
    proportion_covered::PR | false   | _               | _
    dispersal_probs::DP    | false   | _               | _
    human_buffer::B        | false   | _               | _
    dist_buffer::B         | false   | _               | _
    function HumanDispersal{R,W,HP,CS,S,AG,HE,DE,PA,MD,SL,TS,PC,PR,DP,B}(human_pop::HP, cellsize::CS, scale::S, aggregator::AG,
                            human_exponent::HE, dist_exponent::DE, dispersalperpop::PA,
                            max_dispersers::MD, shortlist_len::SL, timestep::TS,
                            dest_shortlists::PC, proportion_covered::PR, dispersal_probs::DP,
                            human_buffer::B, dist_buffer::B
                           ) where {R,W,HP,CS,S,AG,HE,DE,PA,MD,SL,TS,PC,PR,DP,B}

        precalc_human_dispersal!(dest_shortlists, human_pop, cellsize, scale, aggregator,
                                 human_exponent, dist_exponent, shortlist_len, human_buffer, dist_buffer)
        precalc_dispersal_probs!(dispersal_probs, human_pop, dispersalperpop)

        new{R,W,HP,CS,S,AG,HE,DE,PA,MD,SL,TS,PC,PR,DP,B}(human_pop, cellsize, scale, aggregator, human_exponent,
                                 dist_exponent, dispersalperpop, max_dispersers, shortlist_len, timestep, dest_shortlists,
                                 proportion_covered, dispersal_probs, human_buffer, dist_buffer)
    end
end

function HumanDispersal{R,W}(human_pop::HP, cellsize::CS, scale::S, aggregator::AG,
                        human_exponent::HE, dist_exponent::DE, dispersalperpop::PA,
                        max_dispersers::MD, shortlist_len::SL, timestep::TS,
                        dest_shortlists::PC, proportion_covered::PR, dispersal_probs::DP,
                        human_buffer::B, dist_buffer::B
                       ) where {R,W,HP,CS,S,AG,HE,DE,PA,MD,SL,TS,PC,PR,DP,B}
    HumanDispersal{R,W,HP,CS,S,AG,HE,DE,PA,MD,SL,TS,PC,PR,DP,B}(human_pop, cellsize, scale, aggregator, human_exponent,
                             dist_exponent, dispersalperpop, max_dispersers, shortlist_len, timestep, dest_shortlists,
                             proportion_covered, dispersal_probs, human_buffer, dist_buffer)
end

HumanDispersal(; grid=:_default_, human_pop=human_pop, cellsize=1.0, scale=4,
               aggregator=mean, human_exponent=1.0, dist_exponent=1.0, dispersalperpop=1e-3,
               max_dispersers=100.0, shortlist_len=100, timestep=1) = begin

    # Allocate memory
    proportion_covered = nothing
    dispersal_probs = zeros(Union{typeof(dispersalperpop), Missing}, size(human_pop))
    human_buffer = initdownsample(human_pop, scale)
    dist_buffer = initdownsample(human_pop, scale)
    dest_shortlists = init_dest_shortlist(shortlist_len, size(human_buffer))

    HumanDispersal{grid,grid}(human_pop, cellsize, scale, aggregator, human_exponent,
                   dist_exponent, dispersalperpop, max_dispersers, shortlist_len,
                   timestep, dest_shortlists, proportion_covered, dispersal_probs,
                   human_buffer, dist_buffer)
end

getindex(rule::HumanDispersal, I...) = getindex(rule.dest_shortlists, I...)

# Precalculation

"""
Precalculate a dispersal shortlist for each cell
"""
precalc_human_dispersal!(dest_shortlists, human_pop, cellsize, scale, aggregator,
                        human_exponent, dist_exponent, shortlist_len, human_buffer, dist_buffer) = begin
    # Limit shortlist cells to the total available
    @assert shortlist_len <= length(human_buffer)

    # Downsample and get final matrix dimensions
    downsample!(human_buffer, human_pop, aggregator, scale)

    # Precalculate exponentiation of human population matrix
    human_buffer .^= human_exponent

    # Precalculated distances matrix
    for j in 1:size(dist_buffer, 2), i in 1:size(dist_buffer, 1)
        dist_buffer[i, j] = (centercenter_distance(i, j) * cellsize * scale) ^ dist_exponent
    end
    # First cell uses the mean distance from cell centre
    dist_buffer[1, 1] = (cellsize * scale / 6 * (sqrt(2) + log(1 + sqrt(2)))) ^ dist_exponent

    col_mem = [allocate_column(human_buffer, shortlist_len) for thread in 1:Threads.nthreads()]

    Threads.@threads for j in 1:size(human_buffer, 2)
        precalc_col!(dest_shortlists, shortlist_len, human_buffer, dist_buffer, col_mem, j)
    end
end

centercenter_distance(x, y) = sqrt((x - 1)^2 + (y - 1)^2)

init_dest_shortlist(shortlist_len, (h, w)) = init_dest_shortlist(shortlist_len, h, w)
init_dest_shortlist(shortlist_len, h, w) = Precalc[Vector{Interval}(undef, shortlist_len) for i in 1:h, j in 1:w]

allocate_column(human, shortlist_len) = begin
    h, w = size(human)
    gravities = Matrix{Gravity}(undef, size(human)...)
    for j in 1:size(gravities, 2), i in 1:size(gravities, 1)
        gravities[i, j] = Gravity(0.0f0, (i, j))
    end
    gravity_vector = Vector{Gravity}(undef, h * w)
    gravity_shortlist = Vector{Gravity}(undef, shortlist_len)
    interval_shortlist = Vector{Interval}(undef, shortlist_len)
    gravities, gravity_vector, gravity_shortlist, interval_shortlist
end


"""
Precalculate human dispersal shortlist for every cell in a column.
Working on columns is the cleanest way to spread the load accross multiple processors.
"""
function precalc_col!(dest_shortlists, shortlist_len, human, dist, colmem, j)
    h, w = size(human)
    F = eltype(human)
    gravities, gravity_vector, gravity_shortlist, interval_shortlist = colmem[Threads.threadid()]

    for i = 1:h
        cumprop = 0.0f0

        if ismissing(human[i, j])
            dest_shortlists[i, j] = missing
            # Update shortlist proportion matrix to check coverage of the distribution
            # proportion_covered[i, j] = missing
            continue
        end

        # Calculate the gravity for all cells in the grid
        for ii = 1:h, jj = 1:w
            gravities[ii, jj].gravity = if ismissing(human[ii, jj])
                # Missing values are assigned a zero gravity to miss the shortlist
                zero(gravities[1, 1].gravity)
            else
                (human[i, j] * human[ii, jj]) / (dist[abs(i - ii) + 1, abs(j - jj) + 1])
            end
        end

        # Arrange gravities in a vector for 1 dimensional sorting
        for n = 1:size(gravities, 1) * size(gravities, 2)
            gravity_vector[n] = gravities[n]
        end
        # Sort the top shortlist_len gravities in-place, highest first.
        partialsort!(gravity_vector, shortlist_len, rev=true)
        # Copy sorted gravities to the shortlist
        for i in 1:shortlist_len
            gravity_shortlist[i] = gravity_vector[i]
        end

        # Sum gravities in the shortlist
        shortlist_sum::F = sum(gravity_shortlist)

        # Create a list of intervals from the sorted list of gravities.
        # This will be used to choose to randomly select cells from the
        # distribution of gravities
        for (n, m) = enumerate(reverse(gravity_shortlist))
            # Calculaate proportion of current gravity in the complete shortlist
            gravityprop = m.gravity / shortlist_sum
            # Track cumulative proportion for use with `searchsortedfirst()`
            cumprop += gravityprop
            # In case of float innacuracy, finish the distribution at 1.0
            if n == shortlist_len
                cumprop = oneunit(cumprop)
            end
            interval_shortlist[n] = CellInterval(cumprop, m.index)
        end
        # Update output matrix
        dest_shortlists[i, j] .= interval_shortlist

        # Update shortlist proportion matrix to check coverage of the distribution
        # Sum all gravities
        # total_sum::F = sum(gravities)
        # prop = shortlist_sum / total_sum
        # proportion_covered[i, j] = prop
    end
end


"""
Prealculate dispersal probailities for use in the rule

This used to make sense to remove exp(), but probably should be done
on the fly now.
"""
precalc_dispersal_probs!(dispersal_probs, human_activity, dispersalperpop) = begin
    maximum(skipmissing(human_activity)) * dispersalperpop > oneunit(dispersalperpop) &&
        error("dispersalperpop is too high: more propagules can be sent than populaiton")
    dispersal_probs .= human_activity .* dispersalperpop
end

# DynamicGrids Interface
@inline applyrule!(rule::HumanDispersal{R,W}, data, state, index) where {R,W} = begin
    dispersalprob = rule.dispersal_probs[index...]
    ismissing(dispersalprob) && return
    shortlist = rule.dest_shortlists[downsample_index(index, rule.scale)...]
    ismissing(shortlist) && return

    # Find the expected number of dispersers given population, dispersal prob and timeframe
    isnan(state) && println("state", state, " at time: ", currenttime(data))
    isnan(dispersalprob) && println("dispersalprob", dispersalprob)
    total_dispersers = trunc(Int, state * dispersalprob)
    total_dispersers >= zero(total_dispersers) || return

    # Int max number of dispersers in any single dispersal event
    max_dispersers = trunc(Int, rule.max_dispersers)

    # Simulate (possibly) multiple dispersal events from the cell during the timeframe
    dispersed = zero(state)
    while dispersed < total_dispersers
        # Select a subset of the remaining dispersers for a dispersal event
        dispersers = min(rand(1:max_dispersers), total_dispersers - dispersed)
        # Choose a cell to disperse to from the precalculated human dispersal distribution
        dest_id = min(length(shortlist), searchsortedfirst(shortlist, rand()))
        # Randomise cell destination within upsampled destination cells
        upsample = upsample_index(shortlist[dest_id].index, rule.scale)
        dest_index = upsample .+ (rand(0:rule.scale-1), rand(0:rule.scale-1))
        # Skip dispsal to upsampled dest cells that are masked or out of bounds, and try again
        DynamicGrids.ismasked(data, dest_index...) && continue
        DynamicGrids.isinbounds(dest_index, gridsize(data), overflow(data)) || continue
        # Disperse to the cell.
        data[W][dest_index...] += dispersers
        # Track how many have allready dispersed
        dispersed += dispersers
    end
    data[W][index...] -= dispersed
end


# Utilities - removed for memory/performance improvement. Could be returned
# as an optional process.

"""
    populate!(A::AbstractMatrix, rule::HumanDispersal, [I...])
    populate!(A::AbstractMatrix, cells::AbstractArray, [scale=1])

Populate a matrix with the precalculated destinations from a 
[`HumanDispersal`](@ref) rule - either all of the or some subset
if passed the `I...` indexing arguments. This is useful for plotting 
dispersal destinations, especially when used with GeoData.jl
"""
populate!(A::AbstractMatrix, rule::HumanDispersal) = 
    populate!(A, rule.dest_shortlists, rule)
populate!(A::AbstractMatrix, rule::HumanDispersal, I...) =
    populate!(A::AbstractMatrix, rule[I...], rule.scale)
populate!(A::AbstractMatrix, shortlists::AbstractArray, scale=1) = begin
    for I in CartesianIndices(shortlists)
        if ismissing(shortlists[I])
            #A[upsample_index(Tuple(I), rule.scale)...] = missing
        else
            populate!(A, shortlists[I], scale)
        end
    end
    return A
end
populate!(A::AbstractMatrix, cells::Missing, scale) = missing
populate!(A::AbstractMatrix, cells::AbstractVector{CellInterval}, scale=1) = begin
    lastcumprop = 0.0
    for cell in cells
        I = upsample_index(cell.index, scale)
        val = A[I...]
        if ismissing(val)
            A[I...] = cell.cumprop - lastcumprop
        else
            A[I...] += cell.cumprop - lastcumprop
        end
    end
    return A
end

"""
    populate(cells::AbstractVector, size::Tuple, [scale::Int=1])

Returns an array of size `size` populated from the vector
of positions in `cells` rescaled by `scale`.
"""
populate(cells::AbstractVector, size::Tuple, scale::Int=1) = 
    populate!(zeros(Float64, size), cells, scale)

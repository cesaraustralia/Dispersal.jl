
" CellGravity allows sorting on gravity while keeping records of the original cell coordinate "
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
    # fraction::P
    # gravity::M
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
Human driven dispersal rules
"""
abstract type AbstractHumanDispersal <: AbstractPartialRule end

"""
HumanDispersal Rules human-driven dispersal patterns using population data.

Transport connections between grid cells are calculated using distance and human population,
modified with the `human_exponent` and `dist_exponent` parameters. A shortlist of the most
connected cells is selected for use in the simulation.

The time taken for precalulation will depend on the `scale` argument. Values above 1
will downsample the grid to improve precalulation time and runtime performance. A high
scale value is good for use in a live interface.
$(FIELDDOCTABLE)
"""
@description @limits @flattenable struct HumanDispersal{HP,CS,S,AG,HE,DE,EA,MD,SL,TS,PC,PR,DP,B} <: AbstractHumanDispersal
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
    function HumanDispersal{HP,CS,S,AG,HE,DE,PA,MD,SL,TS,PC,PR,DP,B}(human_pop::HP, cellsize::CS, scale::S, aggregator::AG,
                            human_exponent::HE, dist_exponent::DE, dispersalperpop::PA,
                            max_dispersers::MD, shortlist_len::SL, timestep::TS,
                            dest_shortlists::PC, proportion_covered::PR, dispersal_probs::DP,
                            human_buffer::B, dist_buffer::B
                           ) where {HP,CS,S,AG,HE,DE,PA,MD,SL,TS,PC,PR,DP,B}

        precalc_human_dispersal!(dest_shortlists, human_pop, cellsize, scale, aggregator,
                                 human_exponent, dist_exponent, shortlist_len, human_buffer, dist_buffer)
        precalc_dispersal_probs!(dispersal_probs, human_pop, dispersalperpop)

        new{HP,CS,S,AG,HE,DE,PA,MD,SL,TS,PC,PR,DP,B}(human_pop, cellsize, scale, aggregator, human_exponent,
                                 dist_exponent, dispersalperpop, max_dispersers, shortlist_len, timestep, dest_shortlists,
                                 proportion_covered, dispersal_probs, human_buffer, dist_buffer)
    end
end

function HumanDispersal(human_pop::HP, cellsize::CS, scale::S, aggregator::AG,
                        human_exponent::HE, dist_exponent::DE, dispersalperpop::PA,
                        max_dispersers::MD, shortlist_len::SL, timestep::TS,
                        dest_shortlists::PC, proportion_covered::PR, dispersal_probs::DP,
                        human_buffer::B, dist_buffer::B
                       ) where {HP,CS,S,AG,HE,DE,PA,MD,SL,TS,PC,PR,DP,B}
    HumanDispersal{HP,CS,S,AG,HE,DE,PA,MD,SL,TS,PC,PR,DP,B}(human_pop, cellsize, scale, aggregator, human_exponent,
                             dist_exponent, dispersalperpop, max_dispersers, shortlist_len, timestep, dest_shortlists,
                             proportion_covered, dispersal_probs, human_buffer, dist_buffer)
end

HumanDispersal(human_pop; cellsize=1.0, scale=4, aggregator=mean,
               human_exponent=1.0, dist_exponent=1.0, dispersalperpop=1e-3,
               max_dispersers=100.0, shortlist_len=100, timestep=1) = begin

    # Allocate memory
    proportion_covered = nothing
    dispersal_probs = zeros(Union{typeof(dispersalperpop), Missing}, size(human_pop))
    human_buffer = initdownsample(human_pop, scale)
    dist_buffer = initdownsample(human_pop, scale)
    dest_shortlists = init_dest_shortlist(shortlist_len, size(human_buffer))

    HumanDispersal(human_pop, cellsize, scale, aggregator, human_exponent,
                   dist_exponent, dispersalperpop, max_dispersers, shortlist_len,
                   timestep, dest_shortlists, proportion_covered, dispersal_probs,
                   human_buffer, dist_buffer)
end

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
"""
precalc_dispersal_probs!(dispersal_probs, human_activity, dispersalperpop) = begin
    maximum(skipmissing(human_activity)) * dispersalperpop > oneunit(dispersalperpop) && error("dispersalperpop is too high: more propagules can be sent than ppoulaiton")
    dispersal_probs .= human_activity .* dispersalperpop
end



# DynamicGrids Interface

"""
    applyrule(rule::AbstractHumanDispersal, data, state, index)
Simulates human dispersal, weighting dispersal probability based on human
population in the source cell.
"""
applyrule!(rule::AbstractHumanDispersal, data, state, index) = begin
    dispersalprob = rule.dispersal_probs[index...]
    ismissing(dispersalprob) && return

    shortlist = rule.dest_shortlists[downsample_index(index, rule.scale)...]
    ismissing(shortlist) && return

    # Find the expected number of dispersers given population, dispersal prob and timeframe
    isnan(state) && println("state", state, " at time: ", currenttime(data))
    isnan(dispersalprob) && println("dispersalprob", dispersalprob)
    meandispersers = trunc(Int, state * dispersalprob)
    meandispersers >= zero(meandispersers) || return

    # Convert to an actual number of dispersers for this timestep
    # total_dispersers = pois_rand(meandispersers)
    # Check we don't disperse more than the current population (very unlikely with low dispersal probs)
    # if total_dispersers > state
        # total_dispersers = trunc(typeof(total_dispersers), state)
    # end
    total_dispersers = meandispersers # deterministic

    # Get an integer value for the maximum number of dispersers
    # in any single dispersal event
    max_dispersers = trunc(Int, rule.max_dispersers)

    # Simulate (possibly) multiple dispersal events from the cell during the timeframe
    dispersed = 0
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
        DynamicGrids.isinbounds(dest_index, framesize(data), overflow(data)) || continue
        # Disperse to the cell
        data[dest_index...] += dispersers
        # Track how many have allready dispersed
        isnan(dispersers) && error(string("NaN dispersers", (state, index)))
        dispersed += dispersers
    end
    isnan(dispersed) && error(string("NaN disperserd", (state, index)))
    # Subtract dispersed organisms from current cell population
    data[index...] -= dispersed
end


# Utilities

"""
Populate a matrix from a shortlist of cells from one cell in the precalculated matrix
This lets you view the contents of a cell in an AbstractOutput display.

## Arguments:
`a`: A matrix of the same size the precalculation was performed on
`cells`: A vector of [`CellInterval`](@ref)
"""
# populate!(a::AbstractMatrix, cells::Missing, scale) = a
# populate!(a::AbstractMatrix, cells::AbstractVector{<:CellInterval}, scale) = begin
    # for cell in cells
        # a[upsample_index(cell.index, scale)...] = populate_val(a, cell)
    # end
    # a
# end

# populate_val(a::AbstractMatrix{<:Integer}, cell) = 1
# populate_val(a::AbstractMatrix, cell) = cell.fraction

# populate(cells, sze, scale) = populate!(zeros(Float64, sze), cells, scale)

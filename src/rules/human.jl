
"""
CellGravity allows sorting on gravity while keeping records of the original cell coordinate
"""
mutable struct CellGravity{M,I}
    gravity::M
    index::I
end

"""
CellInterval allows ordering a list by the cumulative proportion of the total gravity,
and plotting based on the fraction of total gravity.
"""
struct CellInterval{P,I}
    cumprop::P
    index::I
end

#= Sorting
isless is used to iteratively sort lists and search in funcitons like
`searchsortedfirst` and `partialsort!`. We define `isless` on `CellGravity.gravity` in
order to sortby margnitudes  using `searchsortedfirst`, while tracking the original index.
=#

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
    HumanDispersal{R,W}(mode, human_pop, cellsize, scale, aggregator, human_exponent,
                        dist_exponent, dispersalperpop, max_dispersers, nshortlisted,
                        dest_shortlists, proportion_covered, human_buffer, distances)
    HumanDispersal{R,W}(; mode=BatchGroups(),
                          human_pop,
                          cellsize=1.0,
                          scale=4,
                          aggregator=mean,
                          human_exponent=1.0,
                          dist_exponent=1.0,
                          dispersalperpop=1e-3,
                          max_dispersers=100.0,
                          nshortlisted=100)

Human-driven dispersal patterns using population data.

Transport connections between grid cells are calculated using distance and human population,
modified with the `human_exponent` and `dist_exponent` parameters. A shortlist of the most
connected cells is selected for use in the simulation.

The time taken for precalulation will depend on the `scale` argument. Values above 1
will downsample the grid to improve precalulation time and runtime performance. A high
scale value is good for use in a live interface.

## Keyword Arguments
"""
struct HumanDispersal{R,W,M,HP,CS,S,AG,HE,DE,DP,MD,SL,PC,B,D} <: ManualRule{R,W}
    "Dispersal mode"
    mode::M
    "An array match the grid size containing human population data."
    human_pop::HP
    "The size of the cell width, assuming they are square"
    cellsize::CS
    scale::S
    "A function that aggregates scaled down cells"
    aggregator::AG
    "Human population exponent"
    human_exponent::HE
    "Distance exponent"
    dist_exponent::DE
    "Scales the number of dispersing individuals by human activity (ie population^human_exponent)"
    dispersalperpop::DP
    "Maximum number of dispersers in a dispersal events"
    max_dispersers::MD
    "Length of dest shortlist"
    nshortlisted::SL
    "Array of destination vectors for each cell. Automatically calculated"
    dest_shortlists::PC
    "Buffer array used in precalculation"
    human_buffer::B
    "Buffer array used in precalculation"
    distances::D
    # This constructor is run for every parameter change
    function HumanDispersal{R,W,M,HP,CS,S,AG,HE,DE,DP,MD,SL,PC,B,D}(
        mode::M, human_pop::HP, cellsize::CS, scale::S, aggregator::AG,
        human_exponent::HE, dist_exponent::DE, dispersalperpop::DP,
        max_dispersers::MD, nshortlisted::SL, dest_shortlists::PC,
        human_buffer::B, distances::B
       ) where {R,W,M,HP,CS,S,AG,HE,DE,DP,MD,SL,PC,B,D}

        precalc_human_dispersal!(dest_shortlists, human_pop, cellsize, scale, aggregator,
                                 human_exponent, dist_exponent, nshortlisted, human_buffer, distances)
        new{R,W,M,HP,CS,S,AG,HE,DE,DP,MD,SL,PC,B,D}(mode, human_pop, cellsize, scale, aggregator, human_exponent,
                                 dist_exponent, dispersalperpop, max_dispersers, nshortlisted, dest_shortlists,
                                 human_buffer, distances)
    end
end
# This constructor is run on initialisation with keyword arguments.
# Nothing it defines is affected by rule parameters.
function HumanDispersal{R,W}(;
    mode=BatchGroups(),
    human_pop,
    cellsize=1.0,
    scale=4,
    aggregator=mean,
    human_exponent=Param(1.0; bounds=(1.0, 3.0)),
    dist_exponent=Param(1.0; bounds=(1.0, 3.0)),
    dispersalperpop=Param(1e-3;  bounds=(0.0, 1e-8)) ,
    max_dispersers=Param(100.0, bounds=(50.0, 10000.0)),
    nshortlisted=100,
) where {R,W}
    # Allocate memory
    human_buffer = initdownsample(human_pop, scale)
    distances = initdownsample(human_pop, scale)
    dest_shortlists = alloc_dest_shortlist(nshortlisted, size(human_buffer))

    HumanDispersal{R,W}(
        mode, human_pop, cellsize, scale, aggregator, human_exponent,
        dist_exponent, dispersalperpop, max_dispersers, nshortlisted,
        dest_shortlists, human_buffer, distances
    )
end

# Minimal Array interface

size(rule::HumanDispersal) = size(parent(rule))
getindex(rule::HumanDispersal, I...) = getindex(parent(rule), I...)
parent(rule::HumanDispersal) = rule.dest_shortlists
mode(rule::HumanDispersal) = rule.mode


# Precalculation ###################################################

# Define types used in the precalculation
const Index = Tuple{Int16,Int16}
const Interval = CellInterval{Float32,Index}
const Precalc = Union{Vector{Interval},Missing}
const Gravity = CellGravity{Float32,Index}
const Prop = Union{Float32,Missing}

# Memory preallocation utilities
function alloc_dest_shortlist(nshortlisted, (h, w))
    Precalc[Vector{Interval}(undef, nshortlisted) for i in 1:h, j in 1:w]
end

function alloc_threads(human_buffer, nshortlisted)
    [alloc_gravities(human_buffer, nshortlisted) for thread in 1:Threads.nthreads()]
end

function alloc_gravities(human, nshortlisted)
    gravities = Matrix{Gravity}(undef, size(human)...)
    for j in 1:size(gravities, 2), i in 1:size(gravities, 1)
        gravities[i, j] = Gravity(0.0f0, (i, j))
    end
    gravity_vector = deepcopy(vec(gravities))
    gravities, gravity_vector
end

# Precalculate a dispersal shortlist for each cell
function precalc_human_dispersal!(
    dest_shortlists, human_pop, cellsize, scale, aggregator,
    human_exponent, dist_exponent, nshortlisted, human_buffer, distances
)
    # Limit shortlist cells to the total available
    @assert nshortlisted <= length(human_buffer)

    # Copy downsampled human_pop matrix to human_buffer
    downsample!(human_buffer, human_pop, aggregator, scale)

    precalc_human_exponent!(human_buffer, human_exponent)
    precalc_distances!(distances, dist_exponent, cellsize, scale)

    # Calculate cell gravities by column in separate threads
    thread_alloc = alloc_threads(human_buffer, nshortlisted)
    Threads.@threads for j in 1:size(human_buffer, 2)
        precalc_col!(dest_shortlists, nshortlisted, human_buffer, distances, thread_alloc, j)
    end
end

# Raise the human population/activity matrix to some exponent
function precalc_human_exponent!(human_pop::Matrix, human_exponent::Number)
    human_pop .^= human_exponent
end

# Generate a matrix of distances
function precalc_distances!(distances::Matrix, dist_exponent, cellsize, scale)
    for j in 1:size(distances, 2), i in 1:size(distances, 1)
        distances[i, j] = (centercenter_distance(i, j) * cellsize * scale) ^ dist_exponent
    end
    # First cell uses the mean distance from cell centre
    distances[1, 1] = (cellsize * scale / 6 * (sqrt(2) + log(1 + sqrt(2)))) ^ dist_exponent
    return distances
end

centercenter_distance(a, b) = sqrt((a - 1)^2 + (b - 1)^2)

# Precalculate human dispersal shortlist for every cell in a column.
# Working on columns is a clean way to spread the load accross multiple processors.
function precalc_col!(dest_shortlists, nshortlisted, human, distances, thread_alloc, j)
    col_alloc = thread_alloc[Threads.threadid()]
    for i = 1:size(human, 1)
        precalc_cell!(dest_shortlists, nshortlisted, human, distances, col_alloc, i, j)
    end
end

function precalc_cell!(dest_shortlists, nshortlisted, human, distances, (gravities, gravity_vector), i, j)
    if ismissing(human[i, j])
        dest_shortlists[i, j] = missing
        return
    end
    # Calculate the gravity for all cells in the grid
    assign_gravities!(gravities, human, distances, i, j)
    gravity_vector .= vec(gravities)
    # Sort the top nshortlisted gravities in-place, highest first.
    partialsort!(gravity_vector, nshortlisted, rev=true)
    # Select the top nshortlisted gravities, low to high
    gravity_shortlist = view(gravity_vector, 1:nshortlisted)
    # Convert shortlies gravities to intervals
    gravity2inverval!(dest_shortlists[i, j], gravity_shortlist)
end

# Calculate the gravity for all cells in the grid
function assign_gravities!(gravities, human, distances, i, j)
    T = typeof(first(gravities).gravity)
    for jj = 1:size(human, 2), ii = 1:size(human, 1)
        gravities[ii, jj].gravity = if ismissing(human[ii, jj])
            # Missing values are assigned a zero gravity to miss the shortlist
            zero(T)
        else
            celldist = abs(i - ii) + 1, abs(j - jj) + 1
            (human[ii, jj] / distances[celldist...])
        end
    end
end

# Convert the gravity of a cell to an interval than can
# be randomly selected with `rand()` and `searchsortedfirst`
function gravity2inverval!(interval_shortlist, gravity_shortlist)
    # Sum gravities in the shortlist
    shortlist_sum = sum(gravity_shortlist)
    nshortlisted = length(gravity_shortlist)

    cumprop = 0.0f0
    # Create a list of intervals from the sorted list of gravities.
    # This will be used to randomly choose cells from the
    # distribution of gravities
    for (n, m) = enumerate(reverse(gravity_shortlist))
        # Calculaate proportion of current gravity in the complete shortlist
        gravityprop = m.gravity / shortlist_sum
        # Track cumulative proportion for use with `searchsortedfirst()`
        # In case of float innacuracy, finish the distribution at 1.0
        cumprop = n == nshortlisted ? oneunit(cumprop) : cumprop + gravityprop
        interval_shortlist[n] = CellInterval(cumprop, m.index)
    end
    return interval_shortlist
end


# DynamicGrids Interface ###################################################

@inline function applyrule!(data, rule::HumanDispersal{R,W}, population, cellindex) where {R,W}
    population == zero(population) && return
    dispersalprob = rule.human_pop[cellindex...] * rule.dispersalperpop
    ismissing(dispersalprob) && return
    shortlist = rule.dest_shortlists[downsample_index(cellindex, rule.scale)...]
    #ismissing(shortlist) && return

    #= This formulation introduces a bias where total_dispersers is
    close to max_disersers, for example, for chunks much less than max chunk
    chunks size will allways equal the max chunk size, so the amount wont be random.
    =#
    dispersed = disperse!(
        data[W], mode(rule), rule, shortlist, dispersalprob, population, cellindex
    )
    add!(data[W], -dispersed, cellindex...)
    return 
end


abstract type TransportMode end

struct BatchGroups <: TransportMode end

Base.@kwdef struct HeirarchicalGroups{S} <: TransportMode
    scalar::S = Param(1e-8; bounds=(1e-6, 1e-9))
end

@inline function disperse!(
    data::WritableGridData, mode::HeirarchicalGroups, rule::HumanDispersal,
    shortlist, dispersalprob, population, cellindex
)
    dispersed = zero(population)
    nevents = rand(Binomial(human_pop, dispersalperpop))
    for i in 1:nevents
        maybedispersing = rand(Poissonn(population * mode.scalar))
        # Just exit the loop if we exceed the existing populaiton
        (maybedispersing + dispersed > population) && break
        dispersed += disperse2dest!(data, rule, shortlist, maybedispersing)
    end
    return dispersed
end
@inline function disperse!(
    data::WritableGridData, mode::BatchGroups, rule::HumanDispersal,
    shortlist, dispersalprob, population, cellindex
)
    # Find the expected number of dispersers given population and dispersal prob
    total_dispersers = trunc(Int, min(population * dispersalprob, population))

    # Maximum number of discrete dispersers in any single dispersal event
    # We never know this number, is this defensible
    max_dispersers = trunc(Int, rule.max_dispersers)

    # Simulate (possibly) multiple dispersal events from the cell during the timeframe
    dispersed = zero(population)
    while dispersed < total_dispersers
        # Select a subset of the remaining dispersers for a dispersal event
        maybedispersing = min(rand(1:max_dispersers), total_dispersers - dispersed)
        dispersed += disperse2dest!(data, rule, shortlist, maybedispersing)
        # Track how many have allready dispersed
    end
    return dispersed
end

#= Maybe write the group to a randomised location in the destination array
unless it falls outside the grid or is masked, in which case we
say the event just never happened.
=#
function disperse2dest!(data::DynamicGrids.WritableGridData, rule, shortlist, maybedispersing)
    dest_id = min(length(shortlist), searchsortedfirst(shortlist, rand()))
    # Randomise cell destination within upsampled destination cells
    upsample = upsample_index(shortlist[dest_id].index, rule.scale)
    dest_index = upsample .+ (rand(0:rule.scale-1), rand(0:rule.scale-1))
    # Skip dispsal to upsampled dest cells that are masked or out of bounds, and try again
    return if !DynamicGrids.ismasked(data, dest_index...) &&
        DynamicGrids.isinbounds(dest_index, data)
        # Disperse to the cell.
        add!(data, maybedispersing, dest_index...)
        maybedispersing
    else
        zero(maybedispersing)
    end
end



# Examing the generated shorlists ######################################33

"""
    populate!(A::AbstractMatrix, rule::HumanDispersal, [I...])
    populate!(A::AbstractMatrix, cells::AbstractArray, [scale=1])

Populate a matrix with the precalculated destinations from a
[`HumanDispersal`](@ref) rule - either all of the or some subset
if passed the `I...` indexing arguments. This is useful for plotting
dispersal destinations, especially when used with GeoData.jl
"""
populate!(A::AbstractMatrix, rule::HumanDispersal) =
    populate!(A, rule.dest_shortlists, rule.scale)
populate!(A::AbstractMatrix, rule::HumanDispersal, I...) =
    populate!(A::AbstractMatrix, rule[I...], rule.scale)
# All shortlists for all cells
function populate!(A::AbstractMatrix, shortlists::AbstractArray, scale::Int=1)
    for I in CartesianIndices(shortlists)
        populate!(A, shortlists[I], scale)
    end
    return A
end
# Single shortlist for one cell
function populate!(A::AbstractMatrix, cells::AbstractVector{<:Union{<:CellInterval,Missing}}, scale::Int=1)
    lastcumprop = 0.0
    for cell in cells
        I = upsample_index(cell.index, scale)
        prop = cell.cumprop - lastcumprop
        if ismissing(A[I...])
            A[I...] = prop
        else
            A[I...] += prop
        end
        lastcumprop = cell.cumprop
    end
    return A
end
populate!(A::AbstractMatrix, cells::Missing, scale::Int=1) = missing

"""
    populate(rule::HumanDispersal, [I...])

Returns an array the size of human population matrix filled
with all destination locataion, or with destinations specific
to the passed-in indices `I`.
"""
populate(rule::HumanDispersal, args...) =
    populate!(zeros(Float32, size(rule.human_pop)), rule, args...)

"""
    populate(cells::AbstractVector, size::Tuple, [scale::Int=1])

Returns an array of size `size` populated from the vector
of positions in `cells` rescaled by `scale`.
"""
populate(cells::AbstractVector, size::Tuple, scale::Int=1) =
    populate!(zeros(Float32, size), cells, scale)

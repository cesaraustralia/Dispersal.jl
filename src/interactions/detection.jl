abstract type DetectionModel end

@columns struct ThresholdDetection{T} <: DetectionModel
    # Field      | Default | Flatten | Limits         | Description
    threshold::T | 100.0   | true    | (1.0, 10000.0) | "Number of individuals required in a cell before detection"
end

@columns struct PropabalisticDetection{R,C} <: DetectionModel
    # Field           | Default     | Flatten | Limits       | Description
    detection_rate::R | 0.1         | false   | (0.0, 1.0)   | "Rate of detection per trap"
    trap_coverage::C  | 3.333333e-6 | false   | (0.0, 100.0) | "Proportion of cell coveraged by a single trap"
end

isdetected(d::ThresholdDetection, data, population, ntraps) = population >= d.threshold
isdetected(d::PropabalisticDetection, data, population, ntraps) = begin
    p = 1 - ((1 - d.detection_rate)^ntraps)^(population * d.trap_coverage)
    rand(Binomial(1, p)) == 1
end

@description @limits @flattenable struct Detection{Keys,S,M,T,D,N,J} <: PartialNeighborhoodInteraction{Keys}
    # Field          | Flatten | Limits       | Description
    sites::S         | false   | _            | "Site matrix. Generated from meantraps and sitemask"
    sitemask::M      | false   | _            | "Boolean mask matrix of areas where sites would actually be placed"
    meantraps::T     | true    | (0.0, 100.0) | "Value of mean traps per cell used in generating random trap coverage"
    detectionmode::D | true    | _            | "Model used to determine detection"
    neighborhood::N  | false   | _            | "Response neighborhood. Any `Neighborhood` from DynamicGrids.jl"
    juristictions::J | false   | _            | "Matrix of cells numbered by region"
end
Detection(; population=:population,
            local_response=:local_response,
            juristiction_response=:juristiction_response,
            sites=nothing,
            sitemask=throw(ArgumentError("Must include a mask array ")),
            meantraps=1,
            detectionmode=ThresholdDetection(),
            neighborhood=RadialNeighborhood{1}(),
            juristictions=throw(ArgumentError("Must include a juristictions matrix"))
          ) = begin
    keys = (population, local_response, juristiction_response)
    Detection{keys}(sites, sitemask, meantraps, detectionmode, neighborhood, juristictions)
end
Detection{Keys}(sites::S, sitemask::M, meantraps::T, detectionmode::D, neighborhood::N, juristictions::J) where {Keys,S,M,T,D,N,J} = begin
    sites = build_sites!(sites, sitemask, meantraps)
    Detection{Keys,typeof(sites),M,T,D,N,J}(sites, sitemask, meantraps, detectionmode, neighborhood, juristictions)
end


build_sites!(sites::Nothing, sitemask, meantraps) =
    build_sites!(similar(sitemask, Int8), sitemask, meantraps)
build_sites!(sites::AbstractArray, sitemask, meantraps) = begin
    fill!(sites, false)
    for i = 1:length(sites)
        if !isnothing(mask) && sitemask[i]
            sites[i] = rand(Poisson(meantraps))
        end
    end
    sites
end

radius(i::Detection, ::Val{:local_response}) = radius(neighborhood(i))

@inline applyinteraction!(interaction::Detection{Key}, data::MultiSimData, 
                          (population, local_response, juristiction_response), index) where Key = begin
    POPULATION, LOCAL_RESPONSE, JURISTICTION_RESPONSE = 1, 2, 3
    # TODO what happens on juristiction boundaries...
    local_response && return # Exit it the cell already has a local response
    ntraps = interaction.sites[index...]
    ntraps == 0 && return # Exit if there no traps in the cell
    isdetected(interaction.detectionmode, data, population, ntraps) || return # Exit if nothing detected

    # Start response
    data[LOCAL_RESPONSE][index...] = true
    mapreduceneighbors(setneighbor!, data[LOCAL_RESPONSE], neighborhood(interaction), 
                       interaction, local_response, index)
    sze = size(interaction.juristictions)
    if !juristiction_response
        juristictionid = interaction.juristictions[index...]
        println("set juristiction $juristictionid at $index")
        for j in 1:sze[2], i in 1:sze[1] 
            if interaction.juristictions[i, j] == juristictionid 
                data[JURISTICTION_RESPONSE][i, j] = true
            end
        end
    end
    return
end

# Set neighborhood cells to `true`
@inline setneighbor!(data, hood, interaction::Detection{Keys}, 
                     state, hood_index, dest_index) where Keys =
    @inbounds return data[dest_index...] = true

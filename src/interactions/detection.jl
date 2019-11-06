abstract type DetectionModel end

@columns struct ThresholdDetection{T} <: DetectionModel
    # Field      | Default | Flatten | Limits         | Description
    threshold::T | 100.0   | true    | (1.0, 10000.0) | "Number of individuals required in a cell before detection"
end

isdetected(d::ThresholdDetection, interaction, population) = population >= d.threshold

struct Detection{Keys,S,D,N} <: PartialNeighborhoodInteraction{Keys}
    sites::S
    detection::D
    neighborhood::N
end
Detection(; population=:population,
            quarantine=:quarantine,
            sites=throw(ArgumentError("Must include an array of quarantine sites")),
            detection=ThresholdDetection(),
            neighborhood=RadialNeighborhood{1}()) = begin
    keys = (population, quarantine)
    Detection{keys}(sites, detection, neighborhood)
end
Detection{Keys}(sites::S, detection::D, neighborhood::N) where {Keys,S,D,N} =
    Detection{Keys,S,D,N}(sites, detection, neighborhood)

radius(i::Detection, ::Val{:quarantine}) = radius(neighborhood(i))


@inline applyinteraction!(interaction::Detection{Key}, data::MultiSimData, 
                          (population, quarantine), index) where Key = begin
    # Exit it the cell is already quarantined, or it has no monitoring site
    (quarantine || !interaction.sites[index...]) && return
    # Start quarantine if a population is detected
    if isdetected(interaction.detection, interaction, population)
        data[2][index...] = true
        mapreduceneighbors(setneighbor!, data, neighborhood(interaction), interaction, quarantine, index)
    end
    return
end

# Set neighborhood cells to `true`
@inline setneighbor!(data, hood, interaction::Detection{Keys}, 
                     state, hood_index, dest_index) where Keys =
    @inbounds return data[2][dest_index...] = true

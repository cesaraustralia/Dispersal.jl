abstract type DetectionModel end

struct ThresholdDetection{T} <: DetectionModel
    threshold::T
end

isdetected(d::ThresholdDetection, interaction, pop) = pop > d.threshold

struct Detection{K,D,N} <: Interaction{K}
    detection::D
    neighborhood::N
end
Detection(; organism=:organism,
            quarantine=:quarantine,
            detection=ThresholdDetection(),
            neighborhood=RadialNeighborhood{1}()) = begin
    keys = (organism=organism, quarantine=quarantine)
    Detection{keys}(detection, neighborhood)
end
Detection{K}(detection::D, neighborhood::N) where {K,D,N} = Detection{K,D,N}(neighborhood)

radius(i::Detection{, ::Val{:quarantine}) = radius(i.neighborhood)


@inline applyinteraction!(interaction::Detection{Key}, data, (pop, quarantine), index) where Key = begin
    quarantine && return 
    if isdetected(interaction.detection, interaction, pop)
        DynamicGrids.mapsetneighbor!(data, hood, rule, (pop, quarantine), index)
    end
    return
end

@inline DynamicGrids.setneighbor!(data, hood, interaction::Detection{Keys}, 
                                  state, hood_index, dest_index) where Keys =
    @inbounds return data[Keys[:quarantine]][dest_index...] = true

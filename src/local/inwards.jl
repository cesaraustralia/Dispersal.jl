
"Extend to modify [`InwardsBinaryDispersal`](@ref)"
abstract type AbstractInwardsDispersal <: AbstractNeighborhoodModel end

"""
Binary dispersal within a [`DispersalKernel`](@ref) or other neighborhoods.
Inwards dispersal calculates dispersal *to* the current cell from cells in the neighborhood.
"""
@Probabilistic @Kernel struct InwardsBinaryDispersal{} <: AbstractInwardsDispersal end

@Kernel struct InwardsPopulationDispersal{} <: AbstractInwardsDispersal end

Cellular.radius(model::AbstractInwardsDispersal) = radius(model.neighborhood)
Cellular.temp_neighborhood(model::AbstractInwardsDispersal) = temp_neighborhood(model.neighborhood)

"""
    rule(model::InwardsBinaryDispersal, state::Integer, index, args...)
Runs rule for of [`InwardsBinaryDispersal`](@ref) dispersal.

The current cell is invaded if there is pressure from surrounding cells and
suitable habitat. Otherwise it keeps its current state.
"""
@inline rule(model::InwardsBinaryDispersal, data, state::Integer, args...) = begin
    # Combine neighborhood cells into a single scalar
    cc = neighbors(model.neighborhood, model, data, state, args...)

    # Set to occupied if enough pressure from neighbors
    pressure(model, data.source, cc, args...) ? oneunit(state) : state
end

"""
    rule(model::InwardsPopulationDispersal, state::AbstractFloat, args...)
Runs rule for of [`InwardsPopulationDispersal`](@ref) dispersal.

The current cell is invaded by surrounding cells.
"""
@inline rule(model::InwardsPopulationDispersal, data, state::AbstractFloat, args...) = begin
    p = neighbors(model.neighborhood, model, data, state, args...)
    typeof(state)(rand(Poisson(p)))
end

"""
    neighbors(hood::DispersalKernel, state, index, t, source, dest, args...)

Returns neighbors for a [`DispersalKernel`](@ref), looping over
the array of dispersal propabilities.
"""
@inline neighbors(hood::AbstractDispersalKernel, model::AbstractNeighborhoodModel, data, state,
                  index, args...) = @inbounds return hood.temp ⋅ hood.kernel

"""
    pressure(model, cc)
Calculates the propagule pressure from the output of a neighborhood.
"""
@inline pressure(model, source, cc, args...) = rand() ^ model.prob_threshold > (one(cc) - cc) / one(cc)

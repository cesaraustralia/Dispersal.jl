
"Extend to modify [`InwardsBinaryDispersal`](@ref)"
abstract type AbstractInwardsDispersal <: AbstractNeighborhoodModel end

"""
Binary present/absent dispersal within a [`DispersalKernel`](@ref). 
Inwards dispersal calculates dispersal *to* the current cell from cells in the neighborhood.

The current cell is invaded if there is pressure from surrounding cells and
suitable habitat. Otherwise it keeps its current state.
$(FIELDDOCTABLE)
"""
@Probabilistic @Kernel struct InwardsBinaryDispersal{} <: AbstractInwardsDispersal end

"""
Disperses to the current cells from the populations of the surrounding cells,
using a dispersal kernel.
$(FIELDDOCTABLE)
"""
@Kernel struct InwardsPopulationDispersal{} <: AbstractInwardsDispersal end


"""
Disperses to the current cells from the populations of the surrounding cells,
using a dispersal kernel. Dispersal amounts are randomised with a Poisonn
distribution.
$(FIELDDOCTABLE)
"""
@Kernel struct PoissonInwardsPopulationDispersal{} <: AbstractInwardsDispersal end


Cellular.radius(model::AbstractInwardsDispersal) = radius(model.neighborhood)

Cellular.temp_neighborhood(model::AbstractInwardsDispersal) = temp_neighborhood(model.neighborhood)


@inline rule(model::InwardsBinaryDispersal, data, state::Integer, args...) = begin
    # Combine neighborhood cells into a single scalar
    cc = neighbors(model.neighborhood, model, data, state, args...)

    # Set to occupied if enough pressure from neighbors
    pressure(model, data.source, cc, args...) ? oneunit(state) : state
end

@inline rule(model::InwardsPopulationDispersal, data, state::AbstractFloat, args...) = 
    neighbors(model.neighborhood, model, data, state, args...)

@inline rule(model::PoissonInwardsPopulationDispersal, data, state::AbstractFloat, args...) = begin
    p = neighbors(model.neighborhood, model, data, state, args...)
    p > zero(p) ? typeof(state)(rand(Poisson(p))) : state
end


@inline neighbors(hood::AbstractDispersalKernel, model::AbstractNeighborhoodModel, data, state,
                  index, args...) = @inbounds return hood.temp ⋅ hood.kernel

@inline pressure(model, source, cc, args...) = rand() ^ model.prob_threshold > (one(cc) - cc) / one(cc)

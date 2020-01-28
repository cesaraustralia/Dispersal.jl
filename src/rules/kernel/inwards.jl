
"""
Inwards neighborhood based dispersal models.
"""
abstract type InwardsDispersal <: NeighborhoodRule end

"""
Binary present/absent dispersal within a [`DispersalKernel`](@ref). 
Inwards dispersal calculates dispersal *to* the current cell from cells in the neighborhood.

The current cell is invaded if there is pressure from surrounding cells and
suitable habitat. Otherwise it keeps its current state.
$(FIELDDOCTABLE)
"""
@Probabilistic @Kernel struct InwardsBinaryDispersal{} <: InwardsDispersal end


@inline applyrule(rule::InwardsBinaryDispersal, data, state::Integer, 
                  index, hoodbuffer) = begin
    # Combine neighborhood cells into a single scalar
    s = sumneighbors(neighborhood(rule), hoodbuffer, state)

    # Set to occupied if enough pressure from neighbors
    rand() ^ rule.prob_threshold > (one(s) - s) / one(s) ? oneunit(state) : state
end

"""
Disperses to the current cells from the populations of the 
surrounding cells, using a dispersal kernel deterministically. 

This will only make sense where cell populations are large. 
If they are small, use [PoissonInwardsPopulationDispersal](@ref) to
randomise dispersal jumps.

$(FIELDDOCTABLE)
"""
@Kernel struct InwardsPopulationDispersal{} <: InwardsDispersal end

@inline applyrule(rule::InwardsPopulationDispersal, data, state::AbstractFloat, 
                  index, hoodbuffer) = 
    applykernel(neighborhood(rule), hoodbuffer)

@Layers @Kernel struct SwitchedInwardsPopulationDispersal{Th} <: InwardsDispersal 
    threshold::Th
end

@inline applyrule(rule::SwitchedInwardsPopulationDispersal, data, state::AbstractFloat, 
                  index, hoodbuffer) =
    if layer(rule, data, index) > rule.threshold
        applykernel(neighborhood(rule), hoodbuffer)
    else
        state
    end

DynamicGrids.precalcrules(rule::SwitchedInwardsPopulationDispersal, data) = begin
    precalclayer(layer(rule), rule, data)
end


"""
Disperses to the current cells from the populations of the surrounding cells,
using a dispersal kernel. Dispersal amounts are randomised with a Poisonn
distribution.
$(FIELDDOCTABLE)
"""
@Kernel struct PoissonInwardsPopulationDispersal{} <: InwardsDispersal end

@inline applyrule(rule::PoissonInwardsPopulationDispersal, data, state::AbstractFloat, 
                  index, hoodbuffer) = begin
    p = sumneighbors(neighborhood(rule), hoodbuffer, state)
    p > zero(p) ? typeof(state)(rand(Poisson(p))) : state
end



"Extend to modify [`InwardsBinaryDispersal`](@ref)"
abstract type AbstractInwardsDispersal <: AbstractNeighborhoodRule end

CellularAutomataBase.radius(rule::AbstractInwardsDispersal) = 
    radius(rule.neighborhood)

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




@inline applyrule(rule::InwardsBinaryDispersal, data, state::Integer, args...) = begin
    # Combine neighborhood cells into a single scalar
    cc = neighbors(rule.neighborhood, rule, data, state, args...)

    # Set to occupied if enough pressure from neighbors
    pressure(rule, data.source, cc, args...) ? oneunit(state) : state
end

@inline applyrule(rule::InwardsPopulationDispersal, data, state::AbstractFloat, args...) = 
    neighbors(rule.neighborhood, rule, data, state, args...)

@inline applyrule(rule::PoissonInwardsPopulationDispersal, data, state::AbstractFloat, args...) = begin
    p = neighbors(rule.neighborhood, rule, data, state, args...)
    p > zero(p) ? typeof(state)(rand(Poisson(p))) : state
end


@inline neighbors(hood::AbstractDispersalKernel, rule::AbstractNeighborhoodRule, data, state,
                  index, args...) = @inbounds return buffer(data) ⋅ hood.kernel

@inline pressure(rule, source, cc, args...) = 
    rand() ^ rule.prob_threshold > (one(cc) - cc) / one(cc)

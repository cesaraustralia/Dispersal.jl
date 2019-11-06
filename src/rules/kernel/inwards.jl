
"""
Inwards neighborhood based dispersal models.
"""
abstract type InwardsDispersal{R} <: NeighborhoodRule{R} end

DynamicGrids.radius(rule::InwardsDispersal{R}) where R = R 

# Get the radius from the kernel for all InwardsDispersal
(::Type{T})(kernel, args...) where T <: InwardsDispersal = 
    T{radius(kernel),typeof(kernel),typeof.(args)...}(kernel, args...)

"""
Binary present/absent dispersal within a [`DispersalKernel`](@ref). 
Inwards dispersal calculates dispersal *to* the current cell from cells in the neighborhood.

The current cell is invaded if there is pressure from surrounding cells and
suitable habitat. Otherwise it keeps its current state.
$(FIELDDOCTABLE)
"""
@Kernel @Probabilistic struct InwardsBinaryDispersal{R} <: InwardsDispersal{R} end

@inline applyrule(rule::InwardsBinaryDispersal, data, state::Integer, index, buf) = begin
    # Combine neighborhood cells into a single scalar
    cc = neighbors(neighborhood(rule), rule, buf, state)

    # Set to occupied if enough pressure from neighbors
    pressure(rule, cc) ? oneunit(state) : state
end

"""
Disperses to the current cells from the populations of the surrounding cells,
using a dispersal kernel.
$(FIELDDOCTABLE)
"""
@Kernel struct InwardsPopulationDispersal{R} <: InwardsDispersal{R} end

@inline applyrule(rule::InwardsPopulationDispersal, data, state::AbstractFloat, index, buf) = 
    neighbors(neighborhood(rule), rule, buf, state)


"""
Disperses to the current cells from the populations of the surrounding cells,
using a dispersal kernel. Dispersal amounts are randomised with a Poisonn
distribution.
$(FIELDDOCTABLE)
"""
@Kernel struct PoissonInwardsPopulationDispersal{R} <: InwardsDispersal{R} end

@inline applyrule(rule::PoissonInwardsPopulationDispersal, data, state::AbstractFloat, index, buf) = begin
    p = neighbors(neighborhood(rule), rule, buf, state)
    p > zero(p) ? typeof(state)(rand(Poisson(p))) : state
end


@inline pressure(rule, cc) = rand() ^ rule.prob_threshold > (one(cc) - cc) / one(cc)


"""
Inwards neighborhood based dispersal models.
"""
abstract type InwardsDispersal{R} <: NeighborhoodRule{R} end

DynamicGrids.radius(rule::InwardsDispersal{R}) where R = R 

# Get the radius from the kernel for all InwardsDispersal models
# This simply avoids having to specify R twice.
# It assumes type paremeters and args are both in sequence and allways parametric.
(::Type{T})(args...) where T <: InwardsDispersal = 
    T{findradius(args),typeof.(args)...}(args...)

findradius(args::Tuple{<:Neighborhood,Vararg}) = radius(args[1])
findradius(args::Tuple) = findradius(Base.tail(args))
findradius(::Tuple{}) = 0

"""
Binary present/absent dispersal within a [`DispersalKernel`](@ref). 
Inwards dispersal calculates dispersal *to* the current cell from cells in the neighborhood.

The current cell is invaded if there is pressure from surrounding cells and
suitable habitat. Otherwise it keeps its current state.
$(FIELDDOCTABLE)
"""
@Probabilistic @Kernel struct InwardsBinaryDispersal{R} <: InwardsDispersal{R} end

@inline applyrule(rule::InwardsBinaryDispersal, data, state::Integer, 
                  index, hoodbuffer) = begin
    # Combine neighborhood cells into a single scalar
    cc = neighbors(neighborhood(rule), rule, hoodbuffer, state)

    # Set to occupied if enough pressure from neighbors
    pressure(rule, cc) ? oneunit(state) : state
end

"""
Disperses to the current cells from the populations of the surrounding cells,
using a dispersal kernel.
$(FIELDDOCTABLE)
"""
@Kernel struct InwardsPopulationDispersal{R} <: InwardsDispersal{R} end

@inline applyrule(rule::InwardsPopulationDispersal, data, state::AbstractFloat, 
                  index, hoodbuffer) = 
    applykernel(neighborhood(rule), hoodbuffer)

@Layers @Kernel struct SwitchedInwardsPopulationDispersal{R,Th} <: InwardsDispersal{R} 
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
    rule = precaltimestep(rule, data)
    precalclayer(layer(rule), rule, data)
end


"""
Disperses to the current cells from the populations of the surrounding cells,
using a dispersal kernel. Dispersal amounts are randomised with a Poisonn
distribution.
$(FIELDDOCTABLE)
"""
@Kernel struct PoissonInwardsPopulationDispersal{R} <: InwardsDispersal{R} end

@inline applyrule(rule::PoissonInwardsPopulationDispersal, data, state::AbstractFloat, 
                  index, hoodbuffer) = begin
    p = neighbors(neighborhood(rule), rule, hoodbuffer, state)
    p > zero(p) ? typeof(state)(rand(Poisson(p))) : state
end

@inline pressure(rule, cc) = rand() ^ rule.prob_threshold > (one(cc) - cc) / one(cc)


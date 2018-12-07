"Neighborhoods for dispersal"
abstract type AbstractDispersalNeighborhood <: AbstractNeighborhood end

@chain columns @limits @flattenable @with_kw 


"""
A neighborhood built from a dispersal kernel function and a cell radius.
Can be built directly by passing in the array, radius and overflow
arguments, but preferably use the keyword constructor to build the array from
a dispersal kernel function.
"""
@limits @flattenable struct DispersalNeighborhood{F,P,K,C,I,O} <: AbstractDispersalNeighborhood
    f::F        | false | _
    param::P    | true  | (0.0, 10.0)
    kernel::K   | false | _
    cellsize::C | false | (0.0, 10.0)
    radius::I   | false | (1, 10)
    overflow::O | false | _
    function DispersalNeighborhood{F,P,K,C,I,O}(f::F, param::P, init_kernel::K, 
                                                  cellsize::C, radius::I, overflow::O
                                                 ) where {T,F,P,K,C,I,O}
        kernel = build_dispersal_kernel(f, param, init_kernel, cellsize, radius)
        new{F,P,typeof(kernel),C,I,O}(f, param, kernel, cellsize, radius, overflow)
    end
end

"""
    DispersalNeighborhood(; dir=:inwards, f=exponential, param=1.0, init=[], cellsize=1.0, radius=Int64(3), overflow=Skip())
Constructor for neighborhoods, using a dispersal kernel function and a cell radius.

### Keyword Arguments:
- `f::Function`: any function that accepts a Number argument and returns a propbability between 0.0 and 1.0
- `radius::Integer`: a positive integer
- `overflow = Skip()
"""
DispersalNeighborhood(; f=exponential, param=1.0, init=[], cellsize=1.0, 
                      radius=Int64(3), overflow=Skip()) = begin
    DispersalNeighborhood{typeof.((f, param, init, cellsize, radius, overflow))...
                         }(f, param, init, cellsize, radius, overflow)
end


@mix @columns struct Timestep{TS}
    timestep::TS = 30.0 | true | _
end

@mix struct Neighbors{N}
    "Neighborhood to disperse to or from"
    neighborhood::N = DispersalNeighborhood() | true  | _
end

@mix @columns struct Probabilistic{PT}
    "A real number between one and zero."
    prob_threshold::PT = 0.1 | true | (0.0, 1.0)
end

@mix @columns struct SpotRange{SR}
    "A number or Unitful.jl distance with the same units as cellsize"
    spotrange::SR = 30.0 | true | (0.0, 100.0)
end

@mix @columns struct Fraction{F}
    "The proportion of the population that disperses"
    fraction::F = 0.001 | true | (0.0, 0.1)
end


"Extend to modify [`InwardsBinaryDispersal`](@ref)"
abstract type AbstractInwardsDispersal <: AbstractNeighborhoodModel end

"""
Binary dispersal within a [`DispersalNeighborhood`](@ref) or other neighborhoods.
Inwards dispersal calculates dispersal *to* the current cell from cells in the neighborhood.
"""
@Probabilistic @Neighbors struct InwardsBinaryDispersal{} <: AbstractInwardsDispersal end

@Fraction @Neighbors struct InwardsPopulationDispersal{} <: AbstractInwardsDispersal end


"Extend to modify [`OutwardsBinaryDispersal`](@ref)"
abstract type AbstractOutwardsDispersal <: AbstractPartialNeighborhoodModel end

"""
Binary binary dispersal within a [`DispersalNeighborhood`](@ref)

Outwards dispersal calculates dispersal *from* the current cell to cells
in its neighborhood. This should be more efficient than inwards
dispersal when a small number of cells are occupied, but less efficient when a large
proportion of the grid is occupied.
"""
@Probabilistic @Neighbors struct OutwardsBinaryDispersal{} <: AbstractOutwardsDispersal end

@Fraction @Neighbors struct OutwardsPopulationDispersal{} <: AbstractOutwardsDispersal end


"Extend to modify [`JumpDispersal`](@ref)"
abstract type AbstractJumpDispersal <: AbstractPartialModel end

"Jump dispersal within a [`DispersalNeighborhood`](@ref)] or other neighborhoods."
@Probabilistic @SpotRange struct JumpDispersal{} <: AbstractJumpDispersal end


"""
Abstract supertype that extends `NeighborhoodRule` for neighborhood-based
dispersal rules that update a cell based on the values of the surounding cells,
as if dispersing inwards to the current cell.

The result should be identical to the matching [`OutwardsDispersal`](@ref)
methods.
"""
abstract type InwardsDispersal{R,W} <: NeighborhoodRule{R,W} end

"""
    InwardsBinaryDispersal(neighborhood)
    InwardsBinaryDispersal(; neighborhood=DispersalKernel{3}())
    InwardsBinaryDispersal{R,W}(neighborhood)

Binary present/absent dispersal within a [`DispersalKernel`](@ref).
Inwards dispersal calculates dispersal *to* the current cell from cells in the neighborhood.

The current cell is invaded if there is pressure from surrounding cells and
suitable habitat. Otherwise it keeps its current state.

Pass grid name `Symbol`s to `R` and `W` type parameters to use specific grids.
"""
struct InwardsBinaryDispersal{R,W,NH,PT} <: InwardsDispersal{R,W}
    "Normalised proportions of dispersal to surrounding cells"
    neighborhood::NH
    "A real number between one and zero"
    prob_threshold::PT
end
InwardsBinaryDispersal{R,W}(;
    neighborhood=DispersalKernel{3}(),
    prob_threshold=Param(0.1, bounds=(0.0, 1.0))
) where {R,W} = InwardsBinaryDispersal{R,W}(neighborhood, prob_threshold)


@inline function applyrule(data, rule::InwardsBinaryDispersal, state::Integer, cellindex)
    # Combine neighborhood cells into a single scalar
    s = sum(neighborhood(rule))

    # Set to occupied if enough pressure from neighbors
    rand() ^ rule.prob_threshold > (one(s) - s) / one(s) ? oneunit(state) : state
end

"""
    InwardsPopulationDispersal(neighborhood)
    InwardsPopulationDispersal(; neighborhood=DispersalKernel{3}())
    InwardsPopulationDispersal{R,W}(neighborhood)

Disperses to the current cells from the populations of the
surrounding cells, using a dispersal kernel deterministically.

This will only make sense where cell populations are large.

Pass grid name `Symbol`s to `R` and `W` type parameters to use specific grids.
"""
struct InwardsPopulationDispersal{R,W,NH} <: InwardsDispersal{R,W}
    "Normalised proportions of dispersal to surrounding cells"
    neighborhood::NH
end
InwardsPopulationDispersal{R,W}(; neighborhood=DispersalKernel{3}()) where {R,W} =
    InwardsPopulationDispersal{R,W}(neighborhood)

@inline function applyrule(data, rule::InwardsPopulationDispersal, state, cellindex)
    disperse(neighborhood(rule))
end


struct SwitchedInwardsPopulationDispersal{R,W,NH,A,TI,Th} <: InwardsDispersal{R,W}
    "Normalised proportions of dispersal to surrounding cells"
    neighborhood::NH
    "Key for aux layer"
    auxkey::A
    "Precalculated interpolation indices"
    auxtimeindex::TI
    threshold::Th
end
SwitchedInwardsPopulationDispersal{R,W}(;
    neighborhood=DispersalKernel{3}(),
    auxkey,
    auxtimeindex=1,
    threshold=Param(0.5; bounds=(0.0, 1.0)),
) where {R,W} = SwitchedInwardsPopulationDispersal{R,W}(neighborhood, auxkey, auxtimeindex, threshold)

@inline function applyrule(data, rule::SwitchedInwardsPopulationDispersal, state, cellindex)
    if auxval(data, rule.auxkey, cellindex..., rule.auxtimeindex) > rule.threshold
        disperse(neighborhood(rule))
    else
        state
    end
end

function precalcrule(rule::SwitchedInwardsPopulationDispersal, data)
    precalc_auxtimeindex(aux(data, rule.auxkey), rule, data)
end


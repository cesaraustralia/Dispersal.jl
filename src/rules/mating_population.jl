
struct MatingPopulation{R,W,NH,LC,HC} <: NeighborhoodRule{R,W} 
    "Neighborhood object"
    neighborhood::NH
    "Contribution of local cell to the output"
    local_contribution::LC
    "Contribution of neighbors to the output"
    hood_contribution::HC   
end
function MatingPopulation{R,W}(;
    neighborhood=Moore(1),
    local_contribution=Param(0.5; bounds=(0.0, 1.0)),
    hood_contribution=Param(0.5; bounds=(0.0, 1.0)),
) where {R,W}
    MatingPopulation{R,W}(neighborhood, local_contribution, hood_contribution)
end

# PHENOTYPE and GENOTYPE are the same it's just to update frequency of element

@inline function applyrule(data, rule::MatingPopulation, variable, index)
    @fastmath rule.local_contribution*variable + rule.hood_contribution*mean(neighborhood(rule)) / (rule.local_contribution + rule.hood_contribution)
end

"""
Mating dispersal rules.
"""
struct MatingDispersal{R,W,JX,JY} <: ManualRule{R,W} 
    "Dispersal in X"
    jump_x::JX 
    "Dispersal in Y"
    jump_y::JY 
end
function MatingDispersal{R,W}(;
    jump_x=Param(1; bounds=(1, 100)),
    jump_y=Param(1; bounds=(1, 100)),
) where {R,W}
    MatingDispersal{R,W}(jump_x, jump_y)
end

@inline function applyrule!(data, rule::MatingDispersal{R,W}, state, index) where {R,W}
    # Ignore empty cells
    state > zero(state) || return state
    jump = (rule.jump_x, rule.jump_y)
    jumpdest, is_inbounds = inbounds(jump .+ index, gridsize(data), RemoveOverflow())
    # Update spotted cell if it's on the grid
    if is_inbounds
        @inbounds data[W][jumpdest...] = state
    end
end

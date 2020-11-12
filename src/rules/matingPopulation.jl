@mix @columns struct LocalContribution{LC}
    # Field                   | Default  | Flat  | Bounds     | Description
    localContribution::LC     | 0.5     | true  | (0.0, 1.0)  | "Contribution of local cell to the output"
end

@mix @columns struct HoodContribution{HC}
    # Field                  | Default  | Flat  | Bounds     | Description
    hoodContribution::HC     | 0.5     | true  | (0.0, 1.0)  | "Contribution of neighbors to the output"
end

@Kernel @LocalContribution @HoodContribution struct MatingPopulation{R,W} <: NeighborhoodRule{R,W} end

# PHENOTYPE and GENOTYPE are the same it's just to update frequency of element

@inline applyrule(data, rule::MatingPopulation, variable, index) = begin
    @fastmath rule.localContribution*variable + rule.hoodContribution*mean(neighborhood(rule)) / (rule.localContribution + rule.hoodContribution)
end

"""
Mating dispersal rules.
"""

@mix @columns struct Jump{JX,JY}
    # Field       | Default  | Flat  | Bounds    | Description
    jumpX::JX     | 1        | true  | (1, 100)  | "Dispersal in X"
    jumpY::JY     | 1        | true  | (1, 100)  | "Dispersal in Y"
end

@Jump struct MatingDispersal{R,W} <: ManualRule{R,W} end

@inline applyrule!(data, rule::MatingDispersal{R,W}, state, index) where {R,W} = begin
    # Ignore empty cells
    state > zero(state) || return state
    # 
    jump = (rule.jumpX,rule.jumpY)
    jumpdest, is_inbounds = inbounds(jump .+ index, gridsize(data), RemoveOverflow())
    # Update spotted cell if it's on the grid
    if is_inbounds
        @inbounds data[W][jumpdest...] = state
    end
end
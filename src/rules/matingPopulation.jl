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


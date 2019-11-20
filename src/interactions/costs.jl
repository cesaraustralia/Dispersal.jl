struct Cost{K,Q,TS,TC,I} <: PartialInteraction{K}
    quarantine_cost::Q
    trap_sites::TS
    trap_cost::TC
    industry_value::I
    Cost{K,Q,TS,TC,I}(quarantine_cost::Q, trap_sites::TS, trap_cost::TC, industry_value::I) where {K,Q,TS,TC,I} =
        new{K,Q,TS,TC,I}(quarantine_cost, trap_sites, trap_cost, industry_value)
    Cost{K}(quarantine_cost::Q, trap_sites::TS, trap_cost::TC, industry_value::I) where {K,Q,TS,TC,I} =
        new{K,Q,TS,TC,I}(quarantine_cost, trap_sites, trap_cost, industry_value)
end

Cost(; quarantine=:quarantine,
     cost=:cost,
     population=:population, 
     quarantine_cost=0,
     trap_sites=throw(ArgumentError("Must include an array of trap sites")),
     trap_cost=0,
     industry_value=throw(ArgumentError("Must include an array of annual industry value"))) = begin
    keys = (quarantine, cost, population)
    Cost{keys}(quarantine_cost, trap_sites, trap_cost, industry_value)
end

@inline applyinteraction!(interaction::Cost{Key}, data, (quarantine, cost, population), index) where Key = begin
    QUARANTINE, COST, POPULATION = 1, 2, 3
    if population > 0
        cost += interaction.industry_value[index...] / (10 * 12)
    end
    data[COST][index...] = quarantine * interaction.quarantine_cost + 
                           interaction.trap_sites[index...] * interaction.trap_cost
end

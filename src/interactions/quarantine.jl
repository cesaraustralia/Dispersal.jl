
struct QuarantinedHumanDispersal{K,R,E} <: PartialInteraction{K}
    rule::R
    quarantine_effect::E
    QuarantinedHumanDispersal{K,R,E}(rule::R, quarantine_effect::E) where {K,R,E} = 
        new{K,R,E}(rule, quarantine_effect)
    QuarantinedHumanDispersal{K}(rule::R, quarantine_effect::E) where {K,R,E} = 
        new{K,R,E}(rule, quarantine_effect)
end
QuarantinedHumanDispersal(; population=:population, 
                          quarantine=:quarantine, 
                          rule=HumanDispersal(), 
                          quarantine_effect=1) = begin
    keys = (population, quarantine)
    QuarantinedHumanDispersal{keys}(rule, quarantine_effect)
end


@inline applyinteraction!(interaction::QuarantinedHumanDispersal{Key}, data, state, index) where Key = begin
    POPULATION, QUARANTINE = 1, 2
    rule = interaction.rule
    dispersalprob = rule.dispersal_probs[index...] 
    ismissing(dispersalprob) && return
    quarantined = data[QUARANTINE][index...]
    if quarantined 
        dispersalprob *= interaction.quarantine_effect 
    end
    dispersed = humandispersal!(rule, data[POPULATION], state[POPULATION], index, dispersalprob)
    # These need to stay on separate lines as humandispersal may alter data[index...]
    data[POPULATION][index...] -= dispersed
end

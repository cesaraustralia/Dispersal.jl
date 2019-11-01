struct QuarantinedHumanDispersal{K,D,E} <: Interaction{K}
    dispersalrule::H
    quarantineeffect::E
end
QuarantinedHumanDispersal(; organism=:organism, 
                        quarantine=:quarantine, 
                        dispersalrule=HumanDispersal, 
                        quarantineeffect=1) = begin
    keys = (organism=organism, quarantine=quarantine)
    QuarantinedHumanDispersal{keys}(dispersalrule, management)
end
QuarantinedHumanDispersal{K}(dispersalrule::H, managementeffect::M) where {K,H,M} = 
    QuarantinedHumanDispersal{K,H,M}(dispersalrule, managementeffect)

@inline applyinteraction!(interaction::QuanantinedHumanDispersal{Key}, data, state, index) where Key = begin
    dispersalprob = interaction.rule.dispersal_probs[index...] 
    ismissing(dispersalprob) && return
    quarantined = data[Key[:quarantine]][index...]
    if quarantined 
        dispersalprob *= interaction.quarantineeffect 
    end
    data[Key[:organism]][index...] -= humandispersal!(rule, data, state, index, dispersalprob)
end

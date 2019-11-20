
@description @limits @flattenable struct QuarantinedHumanDispersal{K,R,L,J} <: PartialInteraction{K}
    # Field                | Flatten | Limits     | Description
    rule::R                | false   | _          | "Human dispersal rule"
    local_effect::L        | true    | (0.0, 1.0) | "Scalar 0.0-1.0 modifying the effect of local quarantine on human dispersal"
    juristiction_effect::J | true    | (0.0, 1.0) | "Scalar 0.0-1.0 modifying the effect of juristiction quarantine on human dispersal"
    QuarantinedHumanDispersal{K,R,L,J}(rule::R, local_effect::L, juristiction_effect::J) where {K,R,L,J} = 
        new{K,R,L,J}(rule, local_effect, juristiction_effect)
    QuarantinedHumanDispersal{K}(rule::R, local_effect::L, juristiction_effect::J) where {K,R,L,J} = 
        new{K,R,L,J}(rule, local_effect, juristiction_effect)
end
QuarantinedHumanDispersal(; population=:population, 
                          local_quarantine=:local_quarantine, 
                          juristiction_quarantine=:juristiction_quarantine, 
                          rule=HumanDispersal(), 
                          local_effect=0.1,
                          juristiction_effect=0.1
                         ) = begin
    keys = (population, local_quarantine, juristiction_quarantine)
    QuarantinedHumanDispersal{keys}(rule, local_effect, juristiction_effect)
end

@inline applyinteraction!(interaction::QuarantinedHumanDispersal{Key}, data, 
                          (population, locally_quarantined, juristiction_quarantined), index) where Key = begin
    POPULATION, LOCAL_QUARANTINE, JURISTICTION_QUARANTINE = 1, 2, 3
    rule = interaction.rule
    dispersalprob = rule.dispersal_probs[index...] 
    ismissing(dispersalprob) && return

    shortlist = rule.dest_shortlists[downsample_index(index, rule.scale)...]
    ismissing(shortlist) && return

    if locally_quarantined 
        dispersalprob *= interaction.local_effect 
    end

    # Find the expected number of dispersers given population, dispersal prob and timeframe
    total_dispersers = trunc(Int, population * dispersalprob)
    total_dispersers >= zero(total_dispersers) || return

    # Int max number of dispersers in any single dispersal event
    max_dispersers = trunc(Int, rule.max_dispersers)

    # Simulate (possibly) multiple dispersal events from the cell during the timeframe
    dispersed = zero(population)
    while dispersed < total_dispersers
        # Select a subset of the remaining dispersers for a dispersal event
        dispersers = min(rand(1:max_dispersers), total_dispersers - dispersed)
        # Choose a cell to disperse to from the precalculated human dispersal distribution
        dest_id = min(length(shortlist), searchsortedfirst(shortlist, rand()))
        # Randomise cell destination within upsampled destination cells
        upsample = upsample_index(shortlist[dest_id].index, rule.scale)
        dest_index = upsample .+ (rand(0:rule.scale-1), rand(0:rule.scale-1))
        # Skip dispsal to upsampled dest cells that are masked or out of bounds, and try again
        DynamicGrids.ismasked(data, dest_index...) && continue
        DynamicGrids.isinbounds(dest_index, framesize(data), overflow(data)) || continue
        if juristiction_quarantined && !data[JURISTICTION_QUARANTINE][dest_index...]
            if rand() < interaction.juristiction_effect
                data[POPULATION][dest_index...] += dispersers
            end
        else
            data[POPULATION][dest_index...] += dispersers
        end
        # Track how many have allready dispersed
        dispersed += dispersers
    end
    data[POPULATION][index...] -= dispersed 
end

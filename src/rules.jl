"""
    pressure(model, cc)
Calculates the propagule pressure from the output of a neighborhood.
"""
pressure(model, cc) = rand()^model.prob_threshold > (1 - cc) / 1

"""
    rule(model::AbstractInwardsDispersal, state, index, t, args...)
Runs rule for of [`AbstractInwardsDispersal`](@ref) dispersal.

The current cell is invaded if there is pressure from surrounding cells and
suitable habitat. Otherwise it keeps its current state.
"""
rule(model::AbstractInwardsDispersal, state, index, t, args...) = begin
    # Exit unless cell habitat is suitabile for invasion
    suit = suitability(model.layers, index, t)
    suit >= model.suitability_threshold || return zero(state)

    # Combine neighborhood cells into a single scalar
    cc = neighbors(model.neighborhood, model, state, index, t, args...)

    # Set to occupied if suitable habitat and enough pressure from neighbors
    pressure(model, cc) ? 1 : state
end

"""
    rule(model::AbstractOutwardsDispersal, state, index, t, source, dest, args...)
Runs rule for of [`AbstractOutwardsDispersal`](@ref) dispersal.

Surrounding cells are invaded if the current cell is occupied and they have
suitable habitat. Otherwise they keeps their current state.
"""
rule(model::AbstractOutwardsDispersal, state::Integer, index, t, source, dest, args...) = begin
    state == zero(state) && return # Ignore empty cells 

    propagules = neighbors(model.neighborhood, model, state, index, t, source, dest, args...)

    # Set dest cell state to occupied
    dest[index...] = oneunit(state)
end

rule(model::AbstractOutwardsDispersal, state, index, t, source, dest, args...) = begin
    state == zero(state) && return # Ignore empty cells
    # Grow population - easier to do at the start than the end
    state *= model.growthrate

    propagules = neighbors(model.neighborhood, model, state, index, t, source, dest, args...)

    # Write the new popuation size to the dest array
    dest[row, col] = state - propagules
end

"""
    rule(model::AbstractJumpDispersal, state, index, t, source, dest, args...)
Long range rule for [`AbstractJumpDispersal`](@ref). A random cell
within the spotrange is invaded if it is suitable.
"""
rule(model::AbstractJumpDispersal, state, index, t, source, dest, args...) = begin
    # Ignore empty cells
    state > zero(state) || return
    # Random dispersal events
    rand() < model.prob_threshold || return

    # Calculate maximum spotting distance
    range = -model.spotrange:model.spotrange ./ model.cellsize
    # Randomly select actual spotting distance
    spot = tuple(round.(Int, rand(range, 2) .+ index)...)

    # Update spotted cell if it's on the grid and suitable habitat
    row, col, is_inbounds = inbounds(spot, size(dest), Skip())
    if is_inbounds && suitability(model.layers, (row, col), t) > model.suitability_threshold
        dest[row, col] = oneunit(state)
    end
end

"""
    rule(model::AbstractHumanDispersal, state, index, t, source, dest, args...)
Simulates human dispersal, weighting dispersal probability based on human
population in the source cell.
"""
rule(model::AbstractHumanDispersal, state, index, t, source, dest, args...) = begin
    # Ignore empty cells
    state > zero(state) || return

    rand() < model.prob_threshold * human_impact(model.layers, index, t) || return

    # Calculate maximum spotting distance
    range = -model.spotrange:model.spotrange ./ model.cellsize
    # Randomly select actual spotting distance
    spot = tuple(round.(Int, rand(range, 2) .+ index)...)

    # Update spotted cell if it's on the grid and suitable habitat
    row, col, is_inbounds = inbounds(spot, size(dest), Skip())
    if is_inbounds && suitability(model.layers, (row, col), t) * human_impact(model.layers, (row, col), t) > model.suitability_threshold
        dest[row, col] = oneunit(state)
    end
end

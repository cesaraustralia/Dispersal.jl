"""
    pressure(model, cc)
Calculates the propagule pressure from the output of a neighborhood.
"""
pressure(model, source::Array, cc, args...) = 
    rand() ^ model.prob_threshold > (one(cc) - cc) / one(cc)

"""
    rule(model::AbstractInwardsDispersal, state, row, col, t, layers, args...)
Runs rule for of [`AbstractInwardsDispersal`](@ref) dispersal.

The current cell is invaded if there is pressure from surrounding cells and
suitable habitat. Otherwise it keeps its current state.
"""
rule(model::AbstractInwardsDispersal, state, row, col, t, source, dest, layers, args...) = begin
    # Exit unless cell habitat is suitabile for invasion
    suit = suitability(layers, (row, col), t)
    suit >= model.suitability_threshold || return zero(state)

    # Combine neighborhood cells into a single scalar
    cc = neighbors(model.neighborhood, model, state, row, col, t, source, dest, layers, args...)

    # Set to occupied if suitable habitat and enough pressure from neighbors
    out = pressure(model, source, cc, args...) ? oneunit(state) : state
    out
end

"""
    rule(model::AbstractOutwardsDispersal, state, row, col, t, source, dest, layers, args...)
Runs rule for of [`AbstractOutwardsDispersal`](@ref) dispersal.

Surrounding cells are invaded if the current cell is occupied and they have
suitable habitat. Otherwise they keeps their current state.
"""
rule(model::AbstractOutwardsDispersal, state::Integer, row, col, t, source, dest, layers, args...) = begin
    state == zero(state) && return # Ignore empty cells 

    dest[row, col...] = state 

    propagules = neighbors(model.neighborhood, model, state, row, col, t, source, dest, layers, args...)
end

rule(model::AbstractOutwardsDispersal, state::AbstractFloat, row, col, t, source, dest, layers, args...) = begin
    state == zero(state) && return # Ignore empty cells
    # Grow population - easier to do at the start than the end
    state *= model.growthrate

    propagules = neighbors(model.neighborhood, model, state, row, col, t, source, dest, layers, args...)

    # Write the new popuation size to the dest array
    dest[row, col] = state - propagules
end

"""
    rule(model::AbstractJumpDispersal, state, row, col, t, source, dest, layers, args...)
Long range rule for [`AbstractJumpDispersal`](@ref). A random cell
within the spotrange is invaded if it is suitable.
"""
rule(model::AbstractJumpDispersal, state, row, col, t, source, dest, layers, args...) = begin
    # Ignore empty cells
    state > zero(state) || return
    # Random dispersal events
    spec_rand(source, Float64) < model.prob_threshold || return

    # Calculate maximum spotting distance
    range = -model.spotrange:model.spotrange ./ model.cellsize
    # Randomly select actual spotting distance
    spot = tuple(round.(Int, rand(range, 2) .+ (row, col))...)

    # Update spotted cell if it's on the grid and suitable habitat
    row, col, is_inbounds = inbounds(spot, size(dest), Skip())
    if is_inbounds && suitability(layers, (row, col), t) > model.suitability_threshold
        dest[row, col] = oneunit(state)
    end
end

"""
    rule(model::AbstractHumanDispersal, state, row, col, t, source, dest, args...)
Simulates human dispersal, weighting dispersal probability based on human
population in the source cell.
"""
rule(model::AbstractHumanDispersal, state, row, col, t, source, dest, layers, args...) = begin
    # Ignore empty cells
    state > zero(state) || return

    spec_rand(source, Float64) < model.prob_threshold * human_impact(layers, (row, col), t) || return

    # Calculate maximum spotting distance
    range = -model.spotrange:model.spotrange ./ model.cellsize
    # Randomly select actual spotting distance
    spot = tuple(round.(Int64, rand(range, 2) .+ (row, col))...)

    # Update spotted cell if it's on the grid and suitable habitat
    row, col, is_inbounds = inbounds(spot, size(dest), Skip())
    if is_inbounds#3 && suitability(layers, (row, col), t)# * human_impact(layers, (row, col), t) > model.suitability_threshold
        dest[row, col] = oneunit(state)
    end
end


spec_rand(source::Array, typ, args...) = rand(typ)

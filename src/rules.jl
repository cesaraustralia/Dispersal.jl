"""
    pressure(model, cc)
Calculates the propagule pressure from the output of a neighborhood.
"""
pressure(model, source, cc, args...) = rand() ^ model.prob_threshold > (one(cc) - cc) / one(cc)

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
    pressure(model, source, cc, args...) ? oneunit(state) : state
end

" Grow the population at a fixed rate "
rule(model::FixedRateGrowth, state, args...) = state * model.growthrate


"""
    rule(model::AbstractOutwardsDispersal, state, row, col, t, source, dest, layers, args...)
Runs rule for of [`AbstractOutwardsDispersal`](@ref) dispersal.

Surrounding cells are invaded if the current cell is occupied and they have
suitable habitat. Otherwise they keeps their current state.
"""
rule!(model::AbstractOutwardsDispersal, state::Integer, row, col, t, source, dest, args...) = begin
    state == zero(state) && return # Ignore empty cells 

    propagules = neighbors(model.neighborhood, model, state, row, col, t, source, dest, args...)

    dest[row, col] = state 
    state
end

rule!(model::AbstractOutwardsDispersal, state::AbstractFloat, row, col, t, source, dest, args...) = begin
    state == zero(state) && return state # Ignore empty cells

    propagules = neighbors(model.neighborhood, model, state, row, col, t, source, dest, args...)

    # Write the new popuation size to the dest array
    dest[row, col] = state # - propagules
    state
end

"""
    rule(model::AbstractJumpDispersal, state, row, col, t, source, dest, layers, args...)
Long range rule for [`AbstractJumpDispersal`](@ref). A random cell
within the spotrange is invaded if it is suitable.
"""
rule!(model::AbstractJumpDispersal, state, row, col, t, source, dest, layers, args...) = begin
    # Ignore empty cells
    state > zero(state) || return state
    # Random dispersal events
    spec_rand(source, Float64, args...) < model.prob_threshold || return state

    # Randomly select rpotting distance
    rnge = spec_rand.((source,), (Float64, Float64), tuple.(args)...) .* (model.spotrange / model.cellsize)
    spot = tuple(unsafe_trunc.(Int64, rnge .+ (row, col))...)

    # Update spotted cell if it's on the grid and suitable habitat
    spotrow, spotcol, is_inbounds = inbounds(spot, size(dest), Skip())
    if is_inbounds && suitability(layers, (spotrow, spotcol), t) > model.suitability_threshold
        dest[1, 2] = state # oneunit(state)
    end
    state
end

"""
    rule(model::AbstractHumanDispersal, state, row, col, t, source, dest, args...)
Simulates human dispersal, weighting dispersal probability based on human
population in the source cell.
"""
rule!(model::AbstractHumanDispersal, state, row, col, t, source, dest, layers, args...) = begin
    # Ignore empty cells
    state > zero(state) || return

    spec_rand(source, Float64, args...) < model.prob_threshold * human_impact(layers, (row, col), t) || return

    # Randomly select spotting distance
    rnge = spec_rand.((source,), (Float64, Float64), tuple.(args)...) .* (model.spotrange / model.cellsize)
    spot = tuple(unsafe_trunc.(Int64, rnge .+ (row, col))...)

    # Update spotted cell if it's on the grid and suitable habitat
    spotrow, spotcol, is_inbounds = inbounds(spot, size(dest), Skip())
    if is_inbounds && suitability(layers, (spotrow, spotcol), t) * human_impact(layers, (spotrow, spotcol), t) > model.suitability_threshold
        dest[spotrow, spotcol] = oneunit(state)
    end
    state
end

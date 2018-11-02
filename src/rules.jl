"""
    rule(model::InwardsBinaryDispersal, state::Integer, index, args...)
Runs rule for of [`InwardsBinaryDispersal`](@ref) dispersal.

The current cell is invaded if there is pressure from surrounding cells and
suitable habitat. Otherwise it keeps its current state.
"""
rule(model::InwardsBinaryDispersal, data, state::Integer, index, args...) = begin
    # Combine neighborhood cells into a single scalar
    cc = neighbors(model.neighborhood, model, data, state, index, args...)

    # Set to occupied if enough pressure from neighbors
    pressure(model, data.source, cc, args...) ? oneunit(state) : state
end

"""
    rule(model::InwardsPopulationDispersal, state::AbstractFloat, args...)
Runs rule for of [`InwardsPopulationDispersal`](@ref) dispersal.

The current cell is invaded by surrounding cells.
"""
rule(model::InwardsPopulationDispersal, data, state::AbstractFloat, args...) = 
    state + neighbors(model.neighborhood, model, data, state, args...) * model.fraction


"""
    rule(model::AbstractOutwardsDispersal, data, state, index, layers, args...)
Runs rule for of [`AbstractOutwardsDispersal`](@ref) dispersal.

Surrounding cells are invaded if the current cell is occupied and they have
suitable habitat. Otherwise they keeps their current state.
"""
rule!(model::AbstractOutwardsDispersal, data, state, index, args...) = begin
    state == zero(state) && return # Ignore empty cells 

    propagules = neighbors(model.neighborhood, model, data, state, index, args...)
    
    state
end

"""
    rule(model::AbstractJumpDispersal, data, state, index, layers, args...)
Long range rule for [`AbstractJumpDispersal`](@ref). A random cell
within the spotrange is invaded if it is suitable.
"""
rule!(model::AbstractJumpDispersal, data, state, index, layers, args...) = begin
    # Ignore empty cells
    state > zero(state) || return state

    # Random dispersal events
    spec_rand(data.source, Float64, args...) < model.prob_threshold || return state

    # Randomly select rpotting distance
    rnge = spec_rand.((data.source,), (Float64, Float64), tuple.(args)...) .* (model.spotrange / data.cellsize)
    spot = tuple(unsafe_trunc.(Int64, rnge .+ index)...)
    spot, is_inbounds = inbounds(spot, size(data.dest), Skip())

    # Update spotted cell if it's on the grid
    if is_inbounds
        data.dest[spot...] = state
    end

    state
end

"""
    pressure(model, cc)
Calculates the propagule pressure from the output of a neighborhood.
"""
pressure(model, source, cc, args...) = rand() ^ model.prob_threshold > (one(cc) - cc) / one(cc)


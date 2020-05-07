"""
Abstract supertype. `Objective`s map simulation outputs to predictions 
that can be compared to target data using a loss function.

They must implement [`predictions`](@ref)and [`targets`](@ref) methods.
"""
abstract type Objective end

"""
    targets(obj::Objective)
Returns a targets array given an Objective. The targets must match the size and
dimensions of the prediction array returned by `predictions`.
"""
function targets end

"""
    predictions(obj::Objective, output::Output)

Methods that map an objective object and a simulation output to a
prediction array.
"""
function predictions end

predictions(obj::Objective, output) = output[end]


"""
    SimpleObjective(targets)

A basic objective that holds a target array uses the final frame of the
simulation as the prediction.
"""
struct SimpleObjective{T} <: Objective
    targets::T
end

targets(obj::SimpleObjective) = obj.targets


# """
    # PresenceAbsenceObjective(detectionthreshold, occurance, framesperstep)

# Implementation of a loss objective that converts cell data to regional
# presence/absence and compares to a target of regional occurance data.
# """
# struct PresenceAbsenceObjective{DT,RL,OC,FS,S} <: Objective
#     detectionthreshold::DT
#     occurance::OC
#     framesperstep::FS
# end

# targets(obj::PresenceAbsenceObjective) = obj.occurance
# predictions(obj::PresenceAbsenceObjective, output) = output[end]

"""
Implementation of a loss objective that converts cell data to regional
presence/absence and compares to a target of regional occurance data.

"""
struct RegionObjective{DT,RL,OC,FS,S} <: Objective
    detectionthreshold::DT
    regionlookup::RL
    occurance::OC
    framesperstep::FS
    # TODO fix this offset hack with real DateTime handling
    start::S
end

targets(obj::RegionObjective) = obj.occurance
predictions(obj::RegionObjective, output) = output[end]

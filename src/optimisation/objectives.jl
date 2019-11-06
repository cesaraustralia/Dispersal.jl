"""
Objectives map simulation outputs to predictions that
can be compared to target data using a loss function.

They must implement [`simpredictions`](@ref)and [`targets`](@ref) methods.
"""
abstract type Objective end

"""
    targets(obj::Objective)
Returns a targets array given an Objective. The targets must match the size and
dimensions of the prediction array returned by `simpredictions`.
"""
function targets end

"""
    predictions(obj::Objective, output::Output)
Methods that map an objective object and a simulation output to a
prediction array.
"""
function predictions end



"""
A basic objective that holds a target array uses the final frame of the
simulation as the prediction.
"""
struct SimpleObjective{T} <: Objective
    targets::T
end

targets(obj::SimpleObjective) = obj.targets

predictions(obj::SimpleObjective, output) = output.frames[end]

"""
A simple output that stores each step of the simulation in a vector of arrays.

### Arguments:
- `frames`: Single init array or vector of arrays
- `tstop`: The length of the output.
"""
DynamicGrids.@Output mutable struct RegionOutput{O} <: Output{T}
    objective::O | nothing
end

objective(o::RegionOutput) = o.objective

RegionOutput(objective::Objective; kwargs...) where T = begin
    predictions = [BitArray(zeros(Bool, size(objective.occurance)))]
    RegionOutput(; frames=predictions, objective=objective, kwargs...)
end

DynamicGrids.storeframe!(output::RegionOutput, data::DynamicGrids.SimData, f) = begin
    step = stepfromframe(objective(output), f)
    predictions = output[1]
    for j in 1:framesize(data)[2], i in 1:framesize(data)[1]
        DynamicGrids.blockdo!(data, output, i, j, step, predictions)
    end
end

DynamicGrids.initframes!(output::RegionOutput, init) = begin
    step = stepfromframe(objective(output), 1)
    predictions = output[1]
    predictions .= false
    for j in 1:size(init, 2), i in 1:size(init, 1)
        DynamicGrids.blockdo!(init, output, i, j, step, predictions)
    end
end

@inline DynamicGrids.blockdo!(data, output::RegionOutput, i, j, step, predictions) = begin
    obj = objective(output)
    data[i, j] > obj.detectionthreshold || return
    region = obj.regionlookup[i, j]
    region > zero(region) || return
    predictions[region, step] = true
end

DynamicGrids.showframe(o::RegionOutput, ruleset::Ruleset, f) = nothing

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
predictions(obj::RegionObjective, output::RegionOutput) = output[1]


stepfromframe(objective::RegionObjective, t) = stepfromframe(objective.framesperstep, objective.start, t)
stepfromframe(framesperstep, start, t) = (t - 2one(t) + start) ÷ framesperstep + one(t)

"""
AbstractObjectives map simulation outputs to predictions that
can be compared to target data using a loss function.

They must implement [`simpredictions`](@ref)and [`targets`](@ref) methods.
"""
abstract type AbstractObjective end

"""
    targets(obj::AbstractObjective)
Returns a targets array given an AbstractObjective. The targets must match the size and
dimensions of the prediction array returned by `simpredictions`.
"""
function targets end

"""
    predictions(obj::AbstractObjective, output::AbstractOutput)
Methods that map an objective object and a simulation output to a
prediction array.
"""
function predictions end



"""
A basic objective that holds a target array uses the final frame of the
simulation as the prediction.
"""
struct SimpleObjective{T} <: AbstractObjective
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
DynamicGrids.@Output mutable struct RegionOutput{O} <: AbstractOutput{T}
    objective::O | nothing
end

RegionOutput(frame::AbstractArray{T,2}, starttime, stoptime, objective) where T = begin
    step = stepfromframe(objective, last(DynamicGrids.tspan2fspan(tspan)))
    predictions = [BitArray(zeros(Bool, size(objective.occurance)))]
    RegionOutput{typeof.((predictions, objective))...}(predictions, false, objective)
end

DynamicGrids.storeframe!(output::RegionOutput, data::DynamicGrids.SimData, t) = begin
    step = stepfromframe(output.objective, last(DynamicGrids.tspan2fspan(tspan)))
    predictions = output[1]
    for j in 1:framesize(data)[2], i in 1:framesize(data)[1]
        DynamicGrids.blockdo!(data, output, i, j, step, predictions)
    end
end

DynamicGrids.initframes!(output::RegionOutput, init) = begin
    step = stepfromframe(output.objective, 1)
    predictions = output[1]
    predictions .= false
    for j in 1:size(init, 2), i in 1:size(init, 1)
        DynamicGrids.blockdo!(init, output, i, j, step, predictions)
    end
end

@inline DynamicGrids.blockdo!(data, output::RegionOutput, i, j, step, predictions) = begin
    objective = output.objective
    data[i, j] > objective.detectionthreshold || return
    region = objective.regionlookup[i, j]
    region > zero(region) || return
    predictions[region, step] = true
end

DynamicGrids.showframe(o::RegionOutput, ruleset::AbstractRuleset, t) = nothing

"""
Implementation of a loss objective that converts cell data to regional
presence/absence and compares to a target of regional occurance data.
"""
struct RegionObjective{DT,RL,OC,FS,S} <: AbstractObjective
    detectionthreshold::DT
    regionlookup::RL
    occurance::OC
    framesperstep::FS
    start::S
end

targets(obj::RegionObjective) = obj.occurance
predictions(obj::RegionObjective, output::RegionOutput) = output[1]


stepfromframe(objective::RegionObjective, t) = stepfromframe(objective.framesperstep, objective.start, t)
stepfromframe(framesperstep, start, t) = (t - 2one(t) + start) ÷ framesperstep + one(t)

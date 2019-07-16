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
    simpredictions(obj::AbstractObjective, output::AbstractOutput)
Methods that map an objective object and a simulation output to a 
prediction array.
"""
function simpredictions end



"""
A basic objective that holds a target array uses the final frame of the 
simulation as the prediction.
"""
struct SimpleObjective{T} <: AbstractObjective
    target::T
end

targets(obj::SimpleObjective) = obj.target

simpredictions(obj::SimpleObjective, output) = output.frames[end]



"""
Implementation of a loss objective that converts cell data to regional
presence/absence and compares to a target of regional occurance data.
"""
struct RegionObjective{DT,RL,OC,FS} <: AbstractObjective
    detectionthreshold::DT
    regionlookup::RL
    occurance::OC
    framesperstep::FS
end

targets(obj::RegionObjective) = obj.occurance

simpredictions(obj::RegionObjective, output) = begin
    regions, steps = size(obj.occurance)
    frames = length(output)
    # Allocate arrays for steps and set all cells to zero 
    outputsteps = [similar(output[1]) for f in 1:steps]
    fill!.(outputsteps, zero(eltype(outputsteps[1])))
    # Get the mean population for steps from the frames in each step
    for frame in 1:frames
        step = stepfromframe(obj.framesperstep, frame)
        outputsteps[step] .+= output[frame]
    end
    # Divide all cells by frames per step to get the mean population
    map(s -> s ./= obj.framesperstep, outputsteps)

    # Allocate a boolean array to contain our presence/absence predictions
    prediction = zeros(Bool, size(obj.occurance))
    # Convert mean cell populations to regional prescence/absence
    for t in 1:steps
        for r in 1:regions
            prediction[r, t] = (sum((obj.regionlookup .== r) .& (outputsteps[t] .> 0)) ./
                       sum((obj.regionlookup .== r))) > obj.detectionthreshold
        end
    end
    prediction 
end

stepfromframe(framesperstep, t) = (t - one(t)) ÷ framesperstep + one(t)

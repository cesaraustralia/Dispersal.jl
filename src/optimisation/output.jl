using DynamicGrids: normalise, minval, maxval, ismasked, rgb24

"""
    RegionOutput(init; nframes, objective)

A minimal low-memory output that stores the inhabited regions for
each timestep, as required by the [`RegionObjective`](@ref).
"""
DynamicGrids.@Output mutable struct RegionOutput{O} <: Output{T}
    objective::O | nothing
end
RegionOutput(objective::Objective; kwargs...) where T = begin
    predictions = [BitArray(zeros(Bool, size(objective.occurance)))]
    RegionOutput(; frames=predictions, objective=objective, kwargs...)
end

objective(o::RegionOutput) = o.objective


# RegionObjective and supporting types/methods
#

DynamicGrids.storegrid!(output::RegionOutput, data::DynamicGrids.SimData) = begin
    f = DynamicGrids.gridindex(output, data)
    obj = objective(output)
    step = stepfromframe(obj, f)
    pred = predictions(obj, output)
    for I in CartesianIndices(first(DynamicGrids.init(data)))
        storeprediction!(first(data), output, step, pred, I)
    end
end

DynamicGrids.initgrids!(output::RegionOutput, init) = begin
    obj = objective(output)
    step = stepfromframe(obj, 1)
    pred = predictions(obj, output)
    pred .= false
    for I in CartesianIndices(init)
        storeprediction!(init, output, step, pred, I)
    end
end

"""
Set region presence status in non-zero blocks
"""
@inline storeprediction!(grid, output::RegionOutput, step, predictions, I) = begin
    obj = objective(output)
    grid[I] > obj.detectionthreshold || return
    region = obj.regionlookup[I]
    region > zero(region) || return
    predictions[region, step] = true
    return
end

"""
    ColorRegionFit(objective, truescheme, falsescheme, falsepositivecolor, 
                   truenegativecolor, falsenegativecolor, maskcolor)

An image procesor for visualising the match between predictions and observed 
regional occupancy. A live version of the region fitting optimiser.

## Arguments:

- `objective` : A [`RegionObjective`](@ref) object
- `truescheme` : ColorSchemes.jl scheme or `Greyscale()` for true positive fit, etc.
- `falsescheme` : ColorSchemes.jl scheme or `Greyscale()` for true positive fit, etc.
- `truezerocolor` : `Color`, `Real` beteeen 0.0 and 1.0 or 3 tuple of `Real`.
- `falsezerocolor` : `Color`, `Real` beteeen 0.0 and 1.0 or 3 tuple of `Real`.
- `maskcolor` : `Color` or `Real` between 0.0 and 1.0 to use when a cell is masked.
"""
struct ColorRegionFit{O,T,F,TZ,FZ,M} <: GridProcessor
    objective::O
    truescheme::T
    falsescheme::F
    truezerocolor::TZ
    falsezerocolor::FZ
    maskcolor::M
end

DynamicGrids.grid2image(p::ColorRegionFit, output::ImageOutput, ruleset::Ruleset, grids::NamedTuple, t) = begin
    step = stepfromframe(p.objective, t)
    grid = first(grids)
    img = fill(RGB24(0), size(grid))
    obj = p.objective
    min, max = minval(output), maxval(output)
    for i in CartesianIndices(img)
        region = p.objective.regionlookup[i]
        img[i] = if !(p.maskcolor isa Nothing) && ismasked(mask(ruleset), i) 
            rgb24(p.maskcolor)
        elseif region > zero(region)
            x = grid[i]
            normed = normalise(x, min, max)
            if p.objective.occurance[region, step]
                if !(p.truezerocolor isa Nothing) && normed == zero(normed) 
                    rgb24(p.falsezerocolor)
                else
                    rgb24(p.truescheme, normed)
                end
            else
                if !(p.falsezerocolor isa Nothing) && normed == zero(normed) 
                    rgb24(p.truezerocolor)
                elseif x > obj.detectionthreshold 
                    rgb24(p.falsescheme, normed)
                else
                    rgb24(p.truescheme, normed)
                end
            end
        else
            rgb24(p.maskcolor)
        end
    end
    img
end


stepfromframe(objective::RegionObjective, t) = 
    stepfromframe(objective.framesperstep, objective.start, t)
stepfromframe(framesperstep, start, t) = 
    (t - 2one(t) + start) ÷ framesperstep + one(t)

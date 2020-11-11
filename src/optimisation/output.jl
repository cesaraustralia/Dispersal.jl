using DynamicGrids: normalise, minval, maxval, ismasked, rgb

"""
    RegionOutput(init; tspan, objective, [mask], [extent])

A minimal low-memory output that stores the inhabited regions for
each timestep, as required by the [`RegionObjective`](@ref).
"""
mutable struct RegionOutput{T,F<:AbstractVector{T},E,O} <: Output{T}
    frames::F
    running::Bool
    extent::E
    objective::O
end
RegionOutput(; frames, running, extent, objective, kwargs...) = begin
    predictions = [BitArray(zeros(Bool, size(objective.occurance)))]
    running = false
    RegionOutput(predictions, running, extent, objective)
end

objective(o::RegionOutput) = o.objective


# RegionObjective and supporting types/methods

DynamicGrids.storeframe!(o::RegionOutput, data::DynamicGrids.SimData) = begin
    f = DynamicGrids.frameindex(o, data)
    obj = objective(o)
    step = stepfromframe(obj, f)
    pred = predictions(obj, o)
    for I in CartesianIndices(first(DynamicGrids.init(data)))
        storeprediction!(first(data), o, step, pred, I)
    end
end

DynamicGrids.initgrids!(o::RegionOutput, init::AbstractArray) = begin
    obj = objective(o)
    step = stepfromframe(obj, 1)
    pred = predictions(obj, o)
    pred .= false
    for I in CartesianIndices(init)
        storeprediction!(init, o, step, pred, I)
    end
end
DynamicGrids.initgrids!(o::RegionOutput, init::NamedTuple) =
    DynamicGrids.initgrids!(o, first(init))

"""
Set region presence status in non-zero blocks
"""
@inline storeprediction!(grid, o::RegionOutput, step, predictions, I) = begin
    obj = objective(o)
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

DynamicGrids.grid2image(p::ColorRegionFit, o::ImageOutput, ruleset::Ruleset, grids::NamedTuple, t) = begin
    step = stepfromframe(p.objective, t)
    grid = first(grids)
    img = fill(RGB24(0), size(grid))
    obj = p.objective
    min, max = minval(o), maxval(o)
    for i in CartesianIndices(img)
        region = p.objective.regionlookup[i]
        img[i] = if !(p.maskcolor isa Nothing) && ismasked(mask(o), i)
            rgb(p.maskcolor)
        elseif region > zero(region)
            x = grid[i]
            normed = normalise(x, min, max)
            if p.objective.occurance[region, step]
                if !(p.truezerocolor isa Nothing) && normed == zero(normed)
                    rgb(p.falsezerocolor)
                else
                    rgb(p.truescheme, normed)
                end
            else
                if !(p.falsezerocolor isa Nothing) && normed == zero(normed)
                    rgb(p.truezerocolor)
                elseif x > obj.detectionthreshold
                    rgb(p.falsescheme, normed)
                else
                    rgb(p.truescheme, normed)
                end
            end
        else
            rgb(p.maskcolor)
        end
    end
    img
end


stepfromframe(objective::RegionObjective, t) =
    stepfromframe(objective.framesperstep, objective.start, t)
stepfromframe(framesperstep, start, t) =
    (t - 2one(t) + start) ÷ framesperstep + one(t)

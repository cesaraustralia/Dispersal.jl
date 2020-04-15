using DynamicGrids: normalise, minval, maxval, ismasked, rgb24

"""
    ColorRegionFit(objective, truescheme, falsescheme, falsepositivecolor, 
                   truenegativecolor, falsenegativecolor, maskcolor)

An image procesor for visualising the match between predictions and observed 
regional occupancy. A live version of the region fitting optimiser.

## Arguments:

- `objective` : a RegionObjective object
- `truescheme` : ColorSchemes.jl scheme or `Greyscale()` for true positive fit, etc.
- `falsescheme` : ColorSchemes.jl scheme or `Greyscale()` for true positive fit, etc.
- `truezerocolor` : `Color`, `Real` beteeen 0.0 and 1.0 or 3 tuple of `Real`.
- `falsezerocolor` : `Color`, `Real` beteeen 0.0 and 1.0 or 3 tuple of `Real`.
- `maskcolor` : `Color` or `Real` 0.0 to 1.0 to use when a cell is masked.
"""
struct ColorRegionFit{O,P,N,TZ,FZ,M} <: GridProcessor
    objective::O
    truescheme::P
    falsescheme::N
    truezerocolor::TZ
    falsezerocolor::FZ
    maskcolor::M
end

DynamicGrids.grid2image(p::ColorRegionFit, output::ImageOutput, ruleset::Ruleset, frame, t) = begin
    step = stepfromframe(p.objective, t)
    img = similar(frame, RGB24)
    obj = p.objective
    min, max = minval(ruleset), maxval(ruleset)
    for i in CartesianIndices(frame)
        region = p.objective.regionlookup[i]
        img[i] = if !(p.maskcolor isa Nothing) && ismasked(mask(ruleset), i) 
            p.maskcolor
        elseif region > zero(region)
            x = frame[i]
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
            p.maskcolor
        end
    end
    img
end

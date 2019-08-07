using CellularAutomataBase: normalise, minval, maxval

"""
An image procesor to visualise the model fit, for a live version of
the region fitting optimiser.

Fields:
`objective` : a RegionObjective object
`truepositivecolor` : color of true positive fit, etc.
`falsepositivecolor`
`truenegativecolor`
`falsenegativecolor`
`maskcolor` : color when a cell region of zero or lower
"""
struct ColorRegionFit{O,TP,FP,TN,FN,M} <: AbstractFrameProcessor
    objective::O
    truepositivecolor::TP
    falsepositivecolor::FP
    truenegativecolor::TN
    falsenegativecolor::FN
    maskcolor::M
end

CellularAutomataBase.frametoimage(p::ColorRegionFit, output::AbstractImageOutput, ruleset::AbstractRuleset, frame, t) = begin
    step = stepfromframe(p.objective, t)
    img = similar(frame, RGB24)
    obj = p.objective
    min, max = minval(ruleset), maxval(ruleset)
    for i in CartesianIndices(frame)
        region = p.objective.regionlookup[i]
        img[i] = if region > zero(region)
            x = frame[i]
            if p.objective.occurance[region, step]
                x < obj.detectionthreshold ? rgb(p.falsenegativecolor) : rgb((normalise(x, min, max) .* p.truepositivecolor))
            else
                x < obj.detectionthreshold ? rgb(p.truenegativecolor) : rgb((normalise(x, min, max) .* p.falsepositivecolor))
            end
        else
           rgb(p.maskcolor)
        end
    end
    img
end

rgb(c::RGB24) = c
rgb(c::Tuple) = RGB24(c...)
rgb(c::Number) = RGB24(c)

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

CellularAutomataBase.frametoimage(p::ColorRegionFit, output, frame, t) = begin
    step = stepfromframe(p.objective.framesperstep, t)
    img = similar(frame, RGB24)
    for i in CartesianIndices(frame)
        region = p.objective.regionlookup[i]
        img[i] = if region > zero(region)
            x = frame[i]
            if p.objective.occurance[region, step]
                x == zero(x) ? rgb(p.falsenegativecolor) : rgb((x .* p.truepositivecolor))
            else
                x == zero(x) ? rgb(p.truenegativecolor) : rgb((x .* p.falsepositivecolor))
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

using Cellular: @Ok, @Frames, allocate_frames!, normalize_frame

import Cellular: store_frame!, show_frame, process_image


" An output that condenses a given span of frames to a single frame"
@Ok @Frames struct SumOutput{DI,SF} <: AbstractArrayOutput{T} where DI <: AbstractOutput 
    disp::DI
    frames_per_step::SF
end

SumOutput(frames::AbstractVector, frames_per_step::Number, steps::Number, disp::AbstractOutput) = begin
    output = SumOutput{typeof.((frames, disp, frames_per_step))...}(frames, [false], disp, frames_per_step)
    allocate_frames!(output, frames[1], 2:steps)
    map(f -> f .= 0, output.frames)
    output
end

Cellular.show_frame(output::SumOutput, t::Number) = nothing

" Sums frames on the fly to reduce storage "
Cellular.store_frame!(output::SumOutput, frame, t) = begin
    sze = size(output[1])
    # Determine the timestep being summed to
    ts = step_from_frame(output.frames_per_step, t)
    # Add frame to current sum frame
    for j in 1:sze[2]
        for i in 1:sze[1]
            @inbounds output[ts][i, j] += frame[i, j]
        end
    end
    show_frame(output.disp, output[ts], t)
end

step_from_frame(frames_per_step, t) = (t - one(t)) ÷ frames_per_step + one(t)



" Parametrizer to use with Optim.jl or similar "
struct RegionParametriser{OP,M,I,OC,RL,FS,NR,DT}
    output::OP
    model::M
    init::I
    occurance::OC
    region_lookup::RL
    frames_per_step::FS
    num_replicates::NR
    detection_threshold::DT
end


" Objective function for the parametriser "
(p::RegionParametriser)(params) = begin
    # Rebuild the model with the current parameters
    names = fieldnameflatten(p.model.models)
    println("Parameters: ", collect(zip(names, params)))
    p.model.models = Flatten.reconstruct(p.model.models, params)
    regions, steps = size(p.occurance)
    tstop = steps * p.frames_per_step
    s = zeros(Bool, size(p.occurance))

    cumsum = @distributed (+) for i = 1:p.num_replicates
        o = deepcopy(p.output)
        sim!(o, p.model, p.init; tstop=tstop)
        for t in 1:steps
            for r in 1:regions
                s[r, t] = (sum((p.region_lookup .== r) .& (o[t] .> 0)) ./
                           sum((p.region_lookup .== r))) > p.detection_threshold
            end
        end
        val = sum((s .== p.occurance)) / prod(size(p.occurance))
        println("replicate: ", i, " - accuracy: ", val)
        s
    end
    # then build and array of s array means
    probs = cumsum ./ p.num_replicates
    loss = crossentropy(p.occurance, probs)
    println("cross-entropy loss: ", loss, "\n")

    loss
end

# TODO use LossFunctions.jl so this is modular and other methods can be used
crossentropy(y, p, minprob = 1e-9) = begin
    p = min.(p, 1 - minprob)
    p = max.(p, minprob)
    -sum( y  .* log.(p) .+ (ones(size(y)) .- y) .* log.(ones(size(p)) .- p))
end



""" 
An image procesor to visualise the model fit, for a live version of 
the region fitting optimiser.  

Fields:
`frames_per_step` : The number of frames summed to a single frame
`occurance` : A table of occurrance for each summed step
`region_lookup` : A lookup table matching the occurrance table 
`truepositivecolor` : color of true positive fit, etc.
`falsepositivecolor` 
`truenegativecolor`
`falsenegativecolor`
`maskcolor` : color when a cell region of zero or lower
"""
struct ColorRegionFit{S,OC,CR,TP,FP,TN,FN,M} <: AbstractFrameProcessor 
    frames_per_step::S
    occurance::OC
    region_lookup::CR
    truepositivecolor::TP
    falsepositivecolor::FP
    truenegativecolor::TN
    falsenegativecolor::FN
    maskcolor::M
end

Cellular.process_image(p::ColorRegionFit, output, frame, t) = begin
    step = step_from_frame(p.frames_per_step, t)
    frame = normalize_frame(output, frame)
    img = similar(frame, RGB24) 
    for i in CartesianIndices(frame)
        region = p.region_lookup[i]
        img[i] = if region > zero(region) 
            x = frame[i]
            if p.occurance[region, step]
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

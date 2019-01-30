import Cellular: store_frame!, show_frame, allocate_frames!, 
       @Ok, @Frames, is_async, process_image, normalize_frame


" An output that sums frames for a number of frames "
@Ok @Frames struct SumOutput{DI,SF} <: AbstractArrayOutput{T} where DI <: AbstractOutput disp::DI
    frames_per_step::SF
end

SumOutput(frames::AbstractVector, frames_per_step::Number, steps::Number, disp::AbstractOutput) = begin
    o = SumOutput{typeof.((frames, disp, frames_per_step))...}(frames, [false], disp, frames_per_step)
    allocate_frames!(o, frames[1], 2:steps)
    map(f -> f .= 0, o.frames)
    o
end

# Add Cellular.jl methods

show_frame(output::SumOutput, t::Number) = nothing

" Sums frames on the fly to reduce storage "
store_frame!(o::SumOutput, frame, t) = begin
    sze = size(o[1])
    # Determine the timestep being summed to
    ts = step_from_frame(o.frames_per_step, t)
    # Add frame to current sum frame
    for j in 1:sze[2]
        for i in 1:sze[1]
            @inbounds o[ts][i, j] += frame[i, j]
        end
    end
    show_frame(o.disp, o[ts], t)
end



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

crossentropy(y, p, minprob = 1e-9) = begin
    p = min.(p, 1 - minprob)
    p = max.(p, minprob)
    -sum( y  .* log.(p) .+ (ones(size(y)) .- y) .* log.(ones(size(p)) .- p))
end

step_from_frame(frames_per_step, t) = (t - one(t)) ÷ frames_per_step + one(t)


" An image procesor to visualise the model fit "
struct ColorRegionFit{S,OC,CR,TP,FP,TN,FN,M} <: AbstractImageProcessor 
    frames_per_step::S
    occurance::OC
    region_lookup::CR
    truepositivecolor::TP
    falsepositivecolor::FP
    truenegativecolor::TN
    falsenegativecolor::FN
    maskcolor::M
end

process_image(p::ColorRegionFit, o, frame, t) = begin
    step = step_from_frame(p.frames_per_step, t)
    frame = normalize_frame(o, frame)
    img = similar(frame, RGB24) 
    for i in CartesianIndices(frame)
        region = p.region_lookup[i]
        img[i] = if region > zero(region) 
            x = frame[i]
            if p.occurance[region, step]
                x == zero(x) ? RGB24(p.falsenegativecolor) : RGB24((x .* p.truepositivecolor)...)
            else
                x == zero(x) ? RGB24(p.truenegativecolor) : RGB24((x .* p.falsepositivecolor)...)
            end
        else
           RGB24(p.maskcolor)
        end
    end
    img
end

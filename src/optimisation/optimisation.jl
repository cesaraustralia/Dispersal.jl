using CellularAutomataBase: init

"""
Parametrizer object to use with Optim.jl or similar

# Arguments
`ruleset::Ruleset`: simulation ruleset, with `init` array attached
`objective::AbstractObjective`: objective data
`transform`: single argument function to transform targets and predictions before the loss function
`loss`: LossFunctions.jl loss function
`ngroups`: number of replicate simulation
`groupsize`: number of simulations in a group. Larger groups may inprove distributed performance.
`tstop`: length of simulation
`output`: optional output type. By default an ArrayOutput will be generated.
"""

struct ThreadedReps end
struct DistributedReps end
struct SingleCoreReps end

struct Parametriser{R,OP,OB,TR,L,NR,GS,TS,TH,D,TB,PB}
    ruleset::R
    output::OP
    objective::OB
    transform::TR
    loss::L
    ngroups::NR
    groupsize::GS
    tstop::TS
    threading::TH
    data::D
    targetbuffer::TB
    predictionbuffer::PB
end

Parametriser(ruleset, output, objective, transform, loss, ngroups, groupsize, tstop, threading::ThreadedReps) = begin
    targetbuffer = transform.(targets(objective))
    # Make copies of anything threads will write to
    predictionbuffer = [transform.(targets(objective)) for i in 1:Threads.nthreads()]
    output = [deepcopy(output) for i in 1:Threads.nthreads()]
    data = [CellularAutomataBase.simdata(ruleset, init(ruleset)) for i in 1:Threads.nthreads()]

    Parametriser(ruleset, output, objective, transform, loss, ngroups, groupsize,
                 tstop, threading, data, targetbuffer, predictionbuffer)
end

Parametriser(ruleset, output, objective, transform, loss, ngroups, groupsize, tstop, threading=SingleCoreReps()) = begin
    targetbuffer = transform.(targets(objective))
    predictionbuffer = transform.(targets(objective))
    data = CellularAutomataBase.simdata(ruleset, init(ruleset))

    Parametriser(ruleset, output, objective, transform, loss, ngroups, groupsize,
                 tstop, threading, data, targetbuffer, predictionbuffer)
end

"""
    (p::Parametriser)(params)
Provides an objective function for an optimiser like Optim.jl
"""
(p::Parametriser)(params) = begin
    # Rebuild the rules with the current parameters
    p.ruleset.rules = Flatten.reconstruct(p.ruleset.rules, params, Real)
    obj = p.objective

    names = fieldnameflatten(p.ruleset.rules, Real)
    println("Parameters: \n", ruletypes(typeof(p.ruleset.rules)))
    display(collect(zip(names, params)))

    cumsum = replicate(p.threading, p, params)

    meanloss = cumsum / (p.ngroups * p.groupsize)
    println("mean loss from $(p.ngroups * p.groupsize): ", meanloss, "\n")
    return meanloss
end

replicate(::ThreadedReps, p, params) = begin
    obj = p.objective
    cumsum = Threads.Atomic{Float64}(0.0)
    Threads.@threads for g = 1:p.ngroups
        output = p.output[Threads.threadid()]
        data = p.data[Threads.threadid()]
        predictionbuffer = p.predictionbuffer[Threads.threadid()]
        for n in 1:p.groupsize
            sim!(output, p.ruleset; tstop=p.tstop, data=data)
            loss::Float64 = value(p.loss, p.targetbuffer, predictionbuffer, AggMode.Sum())
            # Core.println("thread: ", Threads.threadid(), " - group: ", g, " - reps ", p.groupsize, " - loss: ", loss)
            Threads.atomic_add!(cumsum, loss)
        end
    end
    cumsum.value
end

replicate(::DistributedReps, p, params) = begin
    obj = p.objective
    cumsum = @distributed (+) for g in 1:p.ngroups
        grouploss = 0.0
        for n in 1:p.groupsize
            sim!(p.output, p.ruleset; tstop=p.tstop)
            p.predictionbuffer .= p.transform.(simpredictions(p.objective, p.output))
            loss::Float64 = value(p.loss, p.targetbuffer, p.predictionbuffer, AggMode.Sum())
            grouploss += loss
        end
        # Core.println("group: ", g, " - reps ($(p.groupsize)) mean loss: ", grouploss/p.groupsize)
        grouploss
    end
    cumsum
end

replicate(::SingleCoreReps, p, params) = begin
    cumsum = 0.0
    obj = p.objective
    for g in 1:p.ngroups
        grouploss = 0.0
        for n in 1:p.groupsize
            sim!(p.output, p.ruleset; tstop=p.tstop, data=p.data)
            p.predictionbuffer .= p.transform.(simpredictions(p.objective, p.output))
            loss::Float64 = value(p.loss, p.targetbuffer, p.predictionbuffer, AggMode.Sum())
            grouploss += loss
        end
        Core.println("group: ", g, " - reps ($(p.groupsize)) mean loss: ", grouploss/p.groupsize)
        cumsum += grouploss
    end
    cumsum
end

struct Accuracy <: LossFunctions.Cost end

LossFunctions.value(::Accuracy, targets, predictions, aggmode) =
    sum(targets .== predictions) / (size(targets, 1) * size(targets, 2))

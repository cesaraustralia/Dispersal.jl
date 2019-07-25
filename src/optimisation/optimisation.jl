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

struct Parametriser{R,OB,F,L,NR,GS,TS,OP,TH}
    ruleset::R
    objective::OB
    transform::F
    loss::L
    ngroups::NR
    groupsize::GS
    tstop::TS
    output::OP
    threading::TH
end

Parametriser(ruleset, objective, transform, loss, ngroups, groupsize, tstop, threading) = begin
    output = ArrayOutput(ruleset.init, tstop)
    Parametriser(ruleset, objective, transform, loss, ngroups, groupsize, tstop, output, threading)
end

"""
    (p::Parametriser)(params)
Provides an objective function for an optimiser like Optim.jl
"""
(p::Parametriser)(params) = begin
    # Rebuild the rules with the current parameters
    p.ruleset.rules = Flatten.reconstruct(p.ruleset.rules, params, Real)

    names = fieldnameflatten(p.ruleset.rules, Real)
    println("Parameters: \n", ruletypes(typeof(p.ruleset.rules)))
    display(collect(zip(names, params)))

    cumsum = replicate(p.threading, p, params)

    meanloss = cumsum / (p.ngroups * p.groupsize)
    println("mean loss from $(p.ngroups * p.groupsize): ", meanloss, "\n")
    return meanloss
end

@inline replicate(::ThreadedReps, p, params) = begin
    cumsum = Threads.Atomic{Float64}(0.0)
    Threads.@threads for g = 1:p.ngroups
        output = p.output[Threads.threadid()]
        for n in 1:p.groupsize
            targs = deepcopy(p.transform.(targets(p.objective)))
            sim!(output, p.ruleset; tstop = p.tstop)
            predictions = p.transform.(simpredictions(p.objective, output))
            loss::Float64 = value(p.loss, targs, predictions, AggMode.Sum())
            Core.println("thread: ", Threads.threadid(), " - group: ", g, " - reps ", p.groupsize, " - loss: ", loss)
            Threads.atomic_add!(cumsum, loss)
        end
    end
    cumsum.value
end

@inline replicate(::DistributedReps, p, params) = begin
    cumsum = @distributed (+) for g in 1:p.ngroups
        grouploss = 0.0
        for n in 1:p.groupsize
            targs = p.transform.(targets(p.objective))
            sim!(p.output, p.ruleset; tstop = p.tstop)
            predictions = p.transform.(simpredictions(p.objective, p.output))
            loss::Float64 = value(p.loss, targs, predictions, AggMode.Sum())
            grouploss += loss
        end
        # Core.println("group: ", g, " - reps ($(p.groupsize)) mean loss: ", grouploss/p.groupsize)
        grouploss
    end
    cumsum
end

@inline replicate(::SingleCoreReps, p, params) = begin
    cumsum = 0.0
    for g in 1:p.ngroups
        grouploss = 0.0
        for n in 1:p.groupsize
            targs = p.transform.(targets(p.objective))
            sim!(p.output, p.ruleset; tstop = p.tstop)
            predictions = p.transform.(simpredictions(p.objective, p.output))
            loss::Float64 = value(p.loss, targs, predictions, AggMode.Sum())
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

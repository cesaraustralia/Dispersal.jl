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
struct Parametriser{R,OB,F,L,NR,GS,TS,OP}
    ruleset::R
    objective::OB
    transform::F
    loss::L
    ngroups::NR
    groupsize::GS
    tstop::TS
    output::OP
end

Parametriser(ruleset, objective, transform, loss, ngroups, groupsize, tstop) = begin
    output = ArrayOutput(ruleset.init, tstop)
    Parametriser(ruleset, objective, transform, loss, ngroups, groupsize, tstop, output)
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

    # cumsum = 0.0
    cumsum = @distributed (+) for i in 1:p.ngroups
        grouploss = 0.0
        for n in 1:p.groupsize
            targs = p.transform.(targets(p.objective))
            sim!(p.output, p.ruleset; tstop = p.tstop)
            predictions = p.transform.(simpredictions(p.objective, p.output))
            loss::Float64 = value(p.loss, targs, predictions, AggMode.Sum())
            grouploss += loss
        end
        Core.println("replicate: ", i, " - group ($(p.groupsize)) mean loss: ", grouploss/p.groupsize)
        # cumsum += grouploss
        grouploss
    end

    meanloss = cumsum / (p.ngroups * p.groupsize)
    println("mean loss from $(p.ngroups * p.groupsize): ", meanloss, "\n")
    return meanloss
end

struct Accuracy <: LossFunctions.Cost end

LossFunctions.value(::Accuracy, targets, predictions, aggmode) =
    sum(targets .== predictions) / (size(targets, 1) * size(targets, 2))


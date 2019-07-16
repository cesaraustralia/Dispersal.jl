""" 
Parametrizer object to use with Optim.jl or similar 

# Arguments
`ruleset::Ruleset`: simulation ruleset, with `init` array attached 
`objective::AbstractObjective`: objective data
`transform`: single argument function to transform targets and predictions before the loss function
`loss`: LossFunctions.jl loss function
`nreplicates`: number of replicate simulation 
`tstop`: length of simulation
`output`: optional output type. By default an ArrayOutput will be generated.
"""
struct Parametriser{R,OB,F,L,NR,TS,OP}
    ruleset::R
    objective::OB
    transform::F
    loss::L
    nreplicates::NR
    tstop::TS
    output::OP
end
Parametriser(ruleset, objective, transform, loss, nreplicates, tstop) = begin
    output = ArrayOutput(ruleset.init, tstop)
    Parametriser(ruleset, objective, transform, loss, nreplicates, tstop, output)
end

""" 
    (p::Parametriser)(params)
Provides an objective function for an optimiser like Optim.jl
"""
(p::Parametriser)(params) = begin
    # Rebuild the rules with the current parameters
    names = fieldnameflatten(p.ruleset.rules, Real)
    println("Parameters: ") 
    println(ruletypes(typeof(p.ruleset.rules)))
    display(collect(zip(names, params)))
    p.ruleset.rules = Flatten.reconstruct(p.ruleset.rules, params, Real)
    i = 1
    targs = p.transform.(targets(p.objective))
    # cumsum = 0.0
    # for i = 1:p.nreplicates
    cumsum = @distributed (+) for i = 1:p.nreplicates
        output = deepcopy(p.output)
        sim!(output, p.ruleset; tstop = p.tstop)
        predictions = p.transform.(simpredictions(p.objective, output))
        loss = value(p.loss, targs, predictions, AggMode.Sum())
        println("replicate: ", i, " - loss: ", loss)
        loss
    end
    meanloss = cumsum / p.nreplicates
    println("mean loss: ", meanloss, "\n")
    return meanloss
end

struct Accuracy <: LossFunctions.Cost end

LossFunctions.value(::Accuracy, targets, predictions, aggmode) = 
    sum(targets .== predictions) / (size(targets, 1) * size(targets, 2))

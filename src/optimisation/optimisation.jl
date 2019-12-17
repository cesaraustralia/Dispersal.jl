using DynamicGrids: init

"Use Thread.@threads to run replicate simulations"
struct ThreadedReplicates end
"Use @distributed to run replicate simulations"
struct DistributedReplicates end
"Use a single processor to run replicate simulations"
struct SingleCoreReplicates end

"""
Parametrizer object to use with Optim.jl or similar

# Arguments
`ruleset::Ruleset`: simulation ruleset, with `init` array attached
`objective::Objective`: objective data
`transform`: single argument function to transform targets and predictions before the loss function
`loss`: LossFunctions.jl loss function
`ngroups`: number of replicate simulation
`groupsize`: number of simulations in a group. Larger groups may inprove distributed performance.
`tstop`: length of simulation
`output`: optional output type. By default an ArrayOutput will be generated.
"""
struct Parametriser{RU,OP,OB,TR,L,NR,GS,TSta,TSto,TH,D,TB,PB,RE}
    ruleset::RU
    output::OP
    objective::OB
    transform::TR
    loss::L
    ngroups::NR
    groupsize::GS
    starttime::TSta
    stoptime::TSto
    threading::TH
    data::D
    targetbuffer::TB
    predictionbuffer::PB
    results::RE
end

DynamicGrids.ruleset(p::Parametriser) = p.ruleset

data(p::Parametriser) = p.data
output(p::Parametriser) = p.output
transform(p::Parametriser) = p.transform
objective(p::Parametriser) = p.objective
loss(p::Parametriser) = p.loss
ngroups(p::Parametriser) = p.ngroups
groupsize(p::Parametriser) = p.groupsize
threading(p::Parametriser) = p.threading
targetbuffer(p::Parametriser) = p.targetbuffer
predictionbuffer(p::Parametriser) = p.predictionbuffer
DynamicGrids.starttime(p::Parametriser) = p.starttime
DynamicGrids.stoptime(p::Parametriser) = p.stoptime
DynamicGrids.tspan(p::Parametriser) = (starttime(p), stoptime(p))


Parametriser(ruleset, output, objective, transform, loss, ngroups, groupsize, starttime, 
             stoptime, threading::ThreadedReplicates) = begin
    targetbuffer = transform.(targets(objective))
    # Make copies of anything threads will write to
    predictionbuffer = [transform.(targets(objective)) for i in 1:Threads.nthreads()]
    output = [deepcopy(output) for i in 1:Threads.nthreads()]
    data = [DynamicGrids.SimData(deepcopy(ruleset), deepcopy(init(ruleset)), starttime) for i in 1:Threads.nthreads()]
    results = [zeros(groupsize) for g in 1:ngroups]

    Parametriser(ruleset, output, objective, transform, loss, ngroups, groupsize,
                 starttime, stoptime, threading, data, targetbuffer, predictionbuffer, results)
end

Parametriser(ruleset, output, objective, transform, loss, ngroups, groupsize, 
             starttime, stoptime, threading=SingleCoreReplicates()) = begin
    targetbuffer = transform.(targets(objective))
    predictionbuffer = transform.(targets(objective))
    data = DynamicGrids.SimData(ruleset, init(ruleset), starttime)
    results = [zeros(groupsize) for g in 1:ngroups]

    Parametriser(ruleset, output, objective, transform, loss, ngroups, groupsize,
                 starttime, stoptime, threading, data, targetbuffer, predictionbuffer, results)
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
    println("\nParameters: \n", ruletypes(typeof(p.ruleset.rules)))
    display(collect(zip(names, params)))
    println()

    replicate!(p.threading, p, params)

    meanloss = sum(sum.(p.results)) / (p.ngroups * p.groupsize)
    println("mean loss from $(p.ngroups * p.groupsize): ", meanloss, "\n")
    return meanloss
end

replicate!(::ThreadedReplicates, p, params) = begin
    obj = p.objective
    Threads.@threads for g in 1:p.ngroups
        id = Threads.threadid()
        output = p.output[id]
        data = p.data[id]
        predictionbuffer = p.predictionbuffer[id]
        for n in 1:p.groupsize
            sim!(output, p.ruleset; tpsan=tspan(p), data=data)
            predictionbuffer .= p.transform.(predictions(obj, output))
            p.results[g][n] = value(p.loss, p.targetbuffer, predictionbuffer, AggMode.Sum())
        end
    end
end

# replicate!(::DistributedReplicates, p, params) = begin
#     obj = p.objective
#     cumsum = @distributed (+) for g in 1:p.ngroups
#         for n in 1:p.groupsize
#             sim!(p.output, p.ruleset; tstop=p.tstop)
#             p.predictionbuffer .= p.transform.(predictions(p.objective, p.output))
#             p.results[g][n] = value(p.loss, p.targetbuffer, p.predictionbuffer, AggMode.Sum())
#         end
#         # Core.println("group: ", g, " - reps ($(p.groupsize)) mean loss: ", grouploss/p.groupsize)
#     end
#     cumsum
# end

replicate!(::SingleCoreReplicates, p, params) = begin
    obj = p.objective
    for g in 1:p.ngroups
        grouploss = 0.0
        for n in 1:p.groupsize
            sim!(p.output, p.ruleset; tspan=tspan(p), data=p.data)
            p.predictionbuffer .= p.transform.(predictions(p.objective, p.output))
            p.results[g][n] = value(p.loss, p.targetbuffer, p.predictionbuffer, AggMode.Sum())
        end
        Core.println("group: ", g, " - reps ($(p.groupsize)) mean loss: ", grouploss/p.groupsize)
    end
end

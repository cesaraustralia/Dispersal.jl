# Rule mixins

@chain columns @default_kw @flattenable @bounds @description 

@mix @columns struct Probabilistic{PT}
    # Field            | Default | Flatn | Limits     | Description
    prob_threshold::PT | 0.1     | true  | (0.0, 1.0) | "A real number between one and zero"
end

@mix @columns struct Timestep{TS,S}
    timestep::TS       | nothing | false | _          | "Timestep converted from sim data. Needs to be separate from rate for DateTime"
    nsteps::S          | 1.0     | false | _          | "The exact nsteps timestep, updated by precalcrule"
end

precalctimestep(rule, data) = precalctimestep(rule.timestep, rule, data)
precalctimestep(ruletimestep::DatePeriod, rule, data) =
    @set rule.nsteps = convert(typeof(rule.nsteps), currenttimestep(data) / Millisecond(ruletimestep))
precalctimestep(ruletimestep::Nothing, rule, data) = @set rule.nsteps = oneunit(rule.nsteps)
precalctimestep(ruletimestep, rule, data) =
    @set rule.nsteps = convert(typeof(rule.nsteps), timestep(data) / ruletimestep)

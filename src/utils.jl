precalc_timestep(rule, data) = precalc_timestep(rule.timestep, rule, data)
precalc_timestep(ruletimestep::DatePeriod, rule, data) =
    @set rule.nsteps = currenttimestep(data) / Millisecond(ruletimestep)
precalc_timestep(ruletimestep::Nothing, rule, data) = @set rule.nsteps = 1
precalc_timestep(ruletimestep, rule, data) =
    @set rule.nsteps = timestep(data) / ruletimestep

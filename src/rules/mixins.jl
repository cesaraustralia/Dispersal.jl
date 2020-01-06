# Rule mixins

@chain columns @description @limits @flattenable @default_kw 

@mix @columns struct Probabilistic{PT}
    # Field            | Default | Flatn | Limits     | Description
    prob_threshold::PT | 0.1     | true  | (0.0, 1.0) | "A real number between one and zero"
end

@mix @columns struct Timestep{TS,S}
    timestep::TS       | nothing | false | _          | "Timestep converted from sim data. Needs to be separate from rate for DateTime"
    nsteps::S          | 1.0     | false | _          | "The exact nsteps timestep, updated by precalcrule"
end

@mix @Timestep struct Layers{L,TI}
    layer::L           | nothing | false | _          | "Data layer"
    timeinterp::TI     | 1       | false | _          | "Precalculated interpolation indices"
end


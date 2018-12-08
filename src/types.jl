"Neighborhoods for dispersal"
abstract type AbstractDispersalNeighborhood <: AbstractNeighborhood end

# Model mixins

@chain columns @limits @flattenable @with_kw 

@mix @columns struct Timestep{TS}
    timestep::TS = 30.0 | true | _
end

@mix @columns struct Probabilistic{PT}
    "A real number between one and zero."
    prob_threshold::PT = 0.1 | true | (0.0, 1.0)
end

@mix @columns struct Fraction{F}
    "The proportion of the population that disperses"
    fraction::F = 0.001 | true | (0.0, 0.1)
end

@mix @columns struct MinMax{M}
    min::M = 0.0       | false | _
    max::M = 1000000.0 | false | _
end

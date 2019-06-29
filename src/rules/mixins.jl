# Rule mixins

@chain columns @description @limits @flattenable @default_kw 

@mix @columns struct Probabilistic{PT}
    # Field            | Def | Flatten | Limits     | Description  
    prob_threshold::PT | 0.1 | true    | (0.0, 1.0) | "A real number between one and zero"
end

@mix @columns struct MinMax{M}
    min::M | 0.0       | false | _ | "Minimum value of the rule"
    max::M | 1000000.0 | false | _ | "Maximum value of the rule"
end

@mix @columns struct Kernel{N}
    neighborhood::N | DispersalKernel() | true  | _ | "Neighborhood to disperse to or from"
end

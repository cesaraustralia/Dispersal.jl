# Rule mixins

@chain columns @description @limits @flattenable @default_kw 

@mix @columns struct Probabilistic{PT}
    # Field            | Def | Flatten | Limits     | Description  
    prob_threshold::PT | 0.1 | true    | (0.0, 1.0) | "A real number between one and zero"
end

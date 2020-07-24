""" 
    AlleeExtinction(minfounders)
    AlleeExtinction(; minfounders=5.0)
    AlleeExtinction{R,W}(minfounders)

Causes extinction in a cell when a population is below the 
minimum number of individuals required to maintain it. 

Pass grid name `Symbol`s to `R` and `W` type parameters to use specific grids.

$(FIELDDOCTABLE)
"""
@columns struct AlleeExtinction{R,W,MF} <: CellRule{R,W}
    # Field         | Default | Flatten | Limits       | Description
    minfounders::MF | 5.0     | true    | (1.0, 200.0) | "Minimum founding individuals required to to start an ongoing population"
end

@inline applyrule(data, model::AlleeExtinction, state, args...) =
    (state >= model.minfounders ? state : zero(state))

""" 
Enforces extinction in a cell with a population below the minimum number of 
individuals required to maintain a population. 
$(FIELDDOCTABLE)
"""
@columns struct AlleeExtinction{MF} <: CellRule
    # Field         | Default | Flatten | Limits       | Description
    minfounders::MF | 5.0     | true    | (1.0, 200.0) | "Minimum founding individuals required to to start an ongoing population"
end

@inline applyrule(model::AlleeExtinction, data, state, args...) =
    (state >= model.minfounders ? state : zero(state))

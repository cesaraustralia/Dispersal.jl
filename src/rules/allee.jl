"Extends AbstractCellRule for allee extinction models"
abstract type AbstractAlleeExtinction <: AbstractCellRule end

""" 
Enforces extinction in a cell without th minimum number of individuals 
required for cell colonisation. 
$(FIELDDOCTABLE)
"""
@columns struct AlleeExtinction{MF} <: AbstractAlleeExtinction
    # Field         | Def | Flatten | Limits       | Description
    minfounders::MF | 5.0 | true    | (1.0, 200.0) | "Minimum founding individuals required to to start an ongoing population"
end

# Rules
@inline applyrule(model::AbstractAlleeExtinction, data, state, args...) =
    (state >= model.minfounders ? state : zero(state))

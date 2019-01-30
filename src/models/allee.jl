# Type declarations
abstract type AbstractAlleeExtinction <: AbstractCellModel end

""" 
Minimum individuals required for cell colonisation 
"""
@columns struct AlleeExtinction{MF} <: AbstractAlleeExtinction
    minfounders::MF = 5.0 | true | (0.0, 200.0)
end

# Rules
@inline rule(model::AbstractAlleeExtinction, data, state, args...) =
    (state >= model.minfounders ? state : zero(state))

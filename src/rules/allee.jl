""" 
    AlleeExtinction(minfounders)
    AlleeExtinction(; minfounders=5.0)
    AlleeExtinction{R,W}(minfounders)

Causes extinction in a cell when a population is below the 
minimum number of individuals required to maintain it. 

Pass grid name `Symbol`s to `R` and `W` type parameters to use specific grids.
"""
struct AlleeExtinction{R,W,MF} <: CellRule{R,W}
    "Minimum founding individuals required to to start an ongoing population"
    minfounders::MF 
end
AlleeExtinction{R,W}(; minfounders=Params(5.0, bounds=(1.0, 200.0))) where {R,W} = 
    AlleeExtinction{R,W}(minfounders)

@inline function applyrule(data, model::AlleeExtinction, state, args...) where {R,W}
    state >= model.minfounders ? state : zero(state)
end

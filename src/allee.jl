""" 
    AlleeExtinction(minfounders)
    AlleeExtinction{R}(; minfounders=5.0)
    AlleeExtinction{R,W}(minfounders)

Causes extinction in a cell when a population is below the 
minimum number of individuals required to maintain it. 

- `minfounders`: minimum founding individuals required to to start an ongoing population.
  Must be a type that can be compared to the grid values using `isless`.

Pass grid `Symbol`s to `R` or both `R` and `W` type parameters to use to specific grids.
"""
struct AlleeExtinction{R,W,MF} <: CellRule{R,W}
    minfounders::MF 
end
AlleeExtinction{R,W}(; minfounders=Param(5.0, bounds=(1.0, 200.0))) where {R,W} = 
    AlleeExtinction{R,W}(minfounders)

@inline function applyrule(data, rule::AlleeExtinction, state, I)
    minfounders = get(data, rule.minfounders, I...)
    state >= minfounders ? state : zero(state)
end

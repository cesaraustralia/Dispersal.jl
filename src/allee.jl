# Type declarations
" Minimum individuals required for cell colonisation "
@columns struct AlleeExtinction{MF} <: AbstractCellModel
    minfounders::MF = 5 | true | (0, 100)
end

# Rules
@inline rule(model::AlleeExtinction, data, state, args...) =
    (state >= model.minfounders ? state : zero(state))

# Type declarations
" Minimum individuals required for cell colonisation "
@columns struct AlleeExtinction{MF} <: AbstractCellModel
    minfounders::MF = 5.0 | true | (0.0, 200.0)
end

# Rules
@inline rule(model::AlleeExtinction, data, state, args...) =
    (state >= model.minfounders ? state : zero(state))

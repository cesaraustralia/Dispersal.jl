# Type declarations
" Minimum individuals required for cell colonisation "
 @columns struct AlleeExtinction{MF} <: AbstractModel
    minfounders::MF = 5 | true | (0.0, 100.0)
end

# Rules
@inline rule(model::AlleeExtinction, data, state, args...) =
    (state > model.minfounders ? state : 0)

# Type declarations

" Minimum individuals required for cell colonisation "
@premix @columns struct MinimumFounders{MF}
    minfounders::MF = 5 | true | (0.0, 100.0)
end

" Simple fixed exponential growth rate solved with Euler method "
@MinimumFounders struct AlleeExtinction{} <: AbstractModel end

# Rules
@inline rule(model::AlleeExtinction, data, state, args...) =
    (state > model.minfounders ? state : 0)

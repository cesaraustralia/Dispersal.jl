using Documenter, Dispersal

makedocs(
    modules = [Dispersal],
    sitename = "Dispersal.jl",
)

deploydocs(
    repo = "github.com/rafaqz/Dispersal.jl.git",
)

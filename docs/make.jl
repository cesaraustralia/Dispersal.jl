using Documenter, DocStringExtensions, Dispersal
using Dispersal: rule, neighbors

makedocs(
    modules = [Dispersal],
    doctest = false,
    clean = false,
    sitename = "Dispersal.jl",
    format = Documenter.HTML(),
    pages = Any[
        "Dispersal" => "index.md",
    ]
)

deploydocs(
    repo = "github.com/rafaqz/Dispersal.jl.git",
    osname = "linux",
    julia = "0.6",
    target = "build",
    deps = nothing,
    make = nothing
)

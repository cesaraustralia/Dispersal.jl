using Documenter, Dispersal, Weave, IJulia, Plots, Pkg, ColorSchemes, Colors, GeoData, JLD2

example = "src/example.jmd"
convert_doc(example, "src/example.md")

makedocs(
    modules = [Dispersal],
    sitename = "Dispersal.jl",
    pages = [
        "Home" => "index.md",
        "Examples" => "example.md"
    ]
)

deploydocs(
    repo = "github.com/rafaqz/Dispersal.jl.git",
)

mkpath.(("build/pdf", "build/notebook"))

weave(example, out_path="build/pdf", doctype="pandoc2pdf")
convert_doc(example, "build/notebook/example.ipynb")

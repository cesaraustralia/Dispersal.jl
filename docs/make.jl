using Documenter, Dispersal, Weave, IJulia

example = "src/example.jmd"
cp(example, "src/example.md", force=true)

# Generate HTML docs
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

pdfdir = "build/pdf" 
notebookdir = "build/notebook"
mkpath.((pdfdir, notebookdir))

# Generate examples pdf
weave(example, out_path=pdfdir, doctype="pandoc2pdf")

# Generate examples notebook
convert_doc(example, joinpath(notebookdir, "example.ipynb"))

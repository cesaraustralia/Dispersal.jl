using Documenter, Dispersal, Weave, IJulia

example = "src/example.jmd"
# Remove YAML
mdlines = readlines(example)
md = join(mdlines[findall(x -> x=="---", mdlines)[2]+1:end])
# Format code blocks for jldoctest 
md = replace(md, Regex("```julia.*") => "```jldoctest")
write("src/example.md", md)

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
    repo = "github.com/cesaraustralia/Dispersal.jl.git",
)

pdfdir = "build/pdf" 
notebookdir = "build/notebook"
mkpath.((pdfdir, notebookdir))

# Generate examples pdf
weave(example, out_path=pdfdir, doctype="pandoc2pdf")

# Generate examples notebook
convert_doc(example, joinpath(notebookdir, "example.ipynb"))

using Documenter, Dispersal, Weave, IJulia

basedir = @__DIR__

example = joinpath(basedir, "src/example.jmd")
# Remove YAML
mdlines = readlines(example)
md = join(mdlines[findall(x -> x=="---", mdlines)[2]+1:end], "\n")
# Format code blocks for jldoctest 
md = replace(md, Regex("```julia.*") => "```jldoctest")
write(joinpath(basedir, "src/example.md"), md)

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

pdfdir = joinpath(basedir, "build/pdf")
notebookdir = joinpath(basedir, "build/notebook")
mkpath.((pdfdir, notebookdir))

# Generate examples pdf
weave(example, out_path=pdfdir, doctype="pandoc2pdf")

# Generate examples notebook
convert_doc(example, joinpath(notebookdir, "example.ipynb"))

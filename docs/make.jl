using Documenter, Dispersal, Weave, IJulia

basedir = @__DIR__

example = joinpath(basedir, "src/example.jmd")

mdpath = joinpath(basedir, "src/example.md")
notebookdir = joinpath(basedir, "build/notebook")

mkpath(joinpath(basedir, "build/assets"))
mkpath(notebookdir)

# Generate examples latex and images
weave(example, out_path=mdpath, doctype="github")

# Generate examples notebook
convert_doc(example, joinpath(notebookdir, "example.ipynb"))

# Generate HTML docs
makedocs(
    modules = [Dispersal],
    sitename = "Dispersal.jl",
    pages = [
        "Home" => "index.md",
        "Examples" => "example.md"
    ],
    clean = false,
)

deploydocs(
    repo = "github.com/cesaraustralia/Dispersal.jl.git",
)

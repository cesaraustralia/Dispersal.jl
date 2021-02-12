using Documenter, Dispersal, Weave, IJulia, Test

basedir = @__DIR__

example = joinpath(basedir, "src/example.jmd")

mdpath = joinpath(basedir, "src/example.md")
notebookdir = joinpath(basedir, "build/notebook")

mkpath(joinpath(basedir, "build/assets"))
mkpath(notebookdir)

# Generate examples latex and images
@test_nowarn weave(example, out_path=mdpath, doctype="github")

# Generate examples notebook
@test_nowarn convert_doc(example, joinpath(notebookdir, "example.ipynb"))

# Generate HTML docs
makedocs(
    modules = [Dispersal],
    sitename = "Dispersal.jl",
    pages = [
        "Home" => "index.md",
        "Examples" => "example.md"
    ],
    clean = false,
    checkdocs = :all,
    strict = true,
)

deploydocs(
    repo = "github.com/cesaraustralia/Dispersal.jl.git",
)

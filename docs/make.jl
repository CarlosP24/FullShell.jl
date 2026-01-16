using Documenter
using FullShell

makedocs(
    sitename = "FullShell.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://CarlosP24.github.io/FullShell.jl",
        assets = String[],
    ),
    modules = [FullShell],
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "API Reference" => "api.md",
    ],
    checkdocs = :exports,
)

deploydocs(
    repo = "github.com/CarlosP24/FullShell.jl.git",
    devbranch = "main",
)

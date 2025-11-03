using ElementaryPDESolutions
using Documenter

DocMeta.setdocmeta!(
    ElementaryPDESolutions, :DocTestSetup,
    :(using ElementaryPDESolutions, StaticArrays);
    recursive = true
)

makedocs(;
    modules = [ElementaryPDESolutions],
    authors = "Thomas G. Anderson, Luiz M. Faria",
    repo = "https://github.com/IntegralEquations/ElementaryPDESolutions.jl/blob/{commit}{path}#{line}",
    sitename = "ElementaryPDESolutions.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://IntegralEquations.github.io/ElementaryPDESolutions.jl",
        edit_link = "main",
        assets = String[]
    ),
    pages = [
        "Home" => "index.md",
        "Docstrings" => "docstrings.md",
    ]
)

deploydocs(;
    repo = "github.com/IntegralEquations/ElementaryPDESolutions.jl",
    devbranch = "main"
)

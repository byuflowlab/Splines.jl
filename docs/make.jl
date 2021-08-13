using Documenter, Splines

makedocs(
    modules = [Splines],
    format = Documenter.HTML(),
    pages = [
        "Intro" => "index.md",
        "Quick Start" => "tutorial.md",
        "Guided Examples" => "howto.md",
        "API Reference" => "reference.md",
        "Theory" => "theory.md",
        "About" => "license.md",
    ],
    repo="https://github.com/byuflowlab/Splines.jl/blob/{commit}{path}#L{line}",
    sitename="Splines.jl",
    authors="Judd Mehr <juddmehr@gmail.com>, Andrew Ning <aning@byu.edu>",
)

deploydocs(
    repo = "github.com/byuflowlab/Splines.jl.git"
)
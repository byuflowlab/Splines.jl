using Documenter, Splines

makedocs(
    modules = [Splines],
    format = Documenter.HTML(),
    pages = ["Home" => "index.md",
             "Functions" => "Functions.md",
             "About" => "license.md"]
    ],
    repo="https://github.com/byuflowlab/Splines.jl/blob/{commit}{path}#L{line}",
    sitename="Splines.jl",
    authors="Judd Mehr <juddmehr@gmail.com>, Andrew Ning <aning@byu.edu>",
)

deploydocs(
    repo = "github.com/byuflowlab/Splies.jl.git"
)
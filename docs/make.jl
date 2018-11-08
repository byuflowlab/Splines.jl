using Documenter, Splines

makedocs(
    root    = ".",
    source  = "./src",
    build   = "./build",
    authors = "Judd Mehr",
    sitename = "Splines.jl",
)

deploydocs(
    repo   = "github.com/byuflowlab/Splines.jl",
    root   = ".",
    target = "./build",
    branch = "gh-pages",
    devbranch = "master",
)
using Documenter, Splines

makedocs(
    root    = ".",
    source  = "./src",
    build   = "./build",
    authors = "Judd Mehr",
    sitename = "Splines.jl",
    pages = ["Home" => "index.md",
             "Functions" => "Functions.md",
             "About" => "license.md"]
)

deploydocs(
    repo   = "github.com/byuflowlab/Splines.jl",
    root   = ".",
    target = "./build",
    branch = "gh-pages",
    devbranch = "master",
)
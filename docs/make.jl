using Documenter, Splines

makedocs(
    root    = ".",
    source  = "./src",
    build   = "./build",
    authors = "Judd Mehr",
    sitename = "Splines.jl",
)

# deploydocs(
#     deps   = Deps.pip("mkdocs", "python-markdown-math"),
#     repo   = "github.com/juddmehr/Splines.jl.git",
#     branch = "gh-pages",
#     latest = "master",
#     julia  = "0.7",
#     osname = "osx"
#     )

deploydocs(
    root   = ".",
    target = "./build",
    repo   = "github.com/juddmehr/Splines.jl.git",
    branch = "gh-pages",
    devbranch = "master",
    osname = "osx"
)
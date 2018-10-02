using Documenter, Splines

makedocs()

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/juddmehr/Splines.jl.git",
    branch = "gh-pages",
    latest = "master",
    julia  = "0.7",
    osname = "osx"
    )
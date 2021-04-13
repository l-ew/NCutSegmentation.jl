using NCutSegmentation
using Documenter

push!(LOAD_PATH, "../src")
makedocs(;
    modules=[NCutSegmentation],
    authors="Luke Ewig <luke.ewig@iob.ch> and contributors",
    repo="https://github.com/l-ew/NCutSegmentation.jl/blob/{commit}{path}#L{line}",
    sitename="NCutSegmentation",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        canonical="https://l-ew.github.io/NCutSegmentation.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    branch="gh-pages",
    devbranch = "main",
    devurl = "stable",
    repo="github.com/l-ew/NCutSegmentation.jl.git",
)

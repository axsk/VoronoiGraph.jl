using VoronoiGraph
using Documenter

DocMeta.setdocmeta!(VoronoiGraph, :DocTestSetup, :(using VoronoiGraph); recursive=true)

makedocs(;
    modules=[VoronoiGraph],
    authors="Sikorski <sikorski@zib.de> and contributors",
    repo="https://github.com/axsk/VoronoiGraph.jl/blob/{commit}{path}#{line}",
    sitename="VoronoiGraph.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://axsk.github.io/VoronoiGraph.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/axsk/VoronoiGraph.jl",
    devbranch="main",
)

push!(LOAD_PATH, "../src/")

using Documenter, ADRIAIndicators

makedocs(
    ;
    sitename="ADRIAIndicators",
    modules = [ADRIAIndicators],
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        assets=["assets/favicon.ico"]
    ),
    pages=[
        "index.md",
        "API" =>
            [
                "metrics.md",
                "cover_metrics.md",
                "juvenile_metrics.md",
                "indices.md",
                "conversions.md"
            ],
    ]
)

deploydocs(;
    repo="github.com/open-AIMS/ADRIAIndicators.jl.git",
    devbranch="main",
    target="build",
    push_preview=false
)

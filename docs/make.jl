push!(LOAD_PATH, "../src/")

using Documenter, ReefMetrics

makedocs(
    ;
    sitename="My Documentation",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        assets=["assets/favicon.ico"]
    ),
    pages=[
        "index.md",
        "API" =>
            [
                "metrics.md",
                "conversions.md",
                "juvenile_metrics.md",
                "indices.md"
            ],
    ]
)

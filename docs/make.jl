push!(LOAD_PATH, "../src/")

using Documenter, ADRIAIndicators

makedocs(
    ;
    sitename="ADRIAIndicators",
    modules = [ADRIAIndicators],
    checkdocs = :public,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true"
    ),
    pages=[
        "Home" => "index.md",
        "usage.md",
        "Cookbook" => [
            "CoralBlox.jl" => "cookbook/coralblox.md",
            "ReefMod Engine" => "cookbook/reefmod_engine.md",
            "C~scape" => "cookbook/cscape.md",
            "CoralPOP" => "cookbook/coralpop.md",
        ],
        "Examples" => [
            "Comparing C~scape and CoralPOP" => "examples/cscape_coralpop.md",
        ],
        "API" =>
            [
                "api/metrics.md",
                "api/cover_metrics.md",
                "api/juvenile_metrics.md",
                "api/indices.md",
                "api/conversions.md"
            ],
        "developing.md",
    ]
)

deploydocs(;
    repo="github.com/open-AIMS/ADRIAIndicators.jl.git",
    devbranch="main",
    target="build",
    push_preview=false
)

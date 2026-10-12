 # Metrics

General metrics that are not specific to a particular group or life stage: indicators of
structural complexity such as `absolute_shelter_volume`, and ecological indicators such as
`coral_diversity`.

```@autodocs
Modules = [ADRIAIndicators]
Pages = joinpath.("src", ["metrics.jl"])
Order   = [:constant, :function, :type]
Private = true
```

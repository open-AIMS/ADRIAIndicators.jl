 # Metrics

This section provides general metrics that are not specific to a particular group or
life stage. These include indicators of structural complexity, such as
`AbsoluteShelterVolume`, and ecological indicators like `CoralDiversity`.

```@autodocs
Modules = [ADRIAIndicators]
Pages = joinpath.("src", ["metrics.jl"])
Order   = [:constant, :function, :type]
Private = true
```

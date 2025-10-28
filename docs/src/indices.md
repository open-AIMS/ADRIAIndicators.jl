# Reef Indices

This section provides high-level reef indices. These indices, such as the
`ReefConditionIndex`, combine multiple lower-level metrics to provide a composite
summary of overall reef state or health for different purposes.

```@autodocs
Modules = [ADRIAIndicators]
Pages = joinpath.("src", ["indices.jl"])
Order   = [:constant, :function, :type]
Private = true
```

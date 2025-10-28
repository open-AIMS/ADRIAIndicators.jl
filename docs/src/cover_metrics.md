# Cover Metrics

This section provides aggregation methods for summarizing coral cover across different
model output dimensions (e.g., time, species, size class). These functions help
reduce high-dimensional coral cover data into simpler, more interpretable summary
metrics.

```@autodocs
Modules = [ADRIAIndicators]
Pages = joinpath.("src", ["cover_metrics.jl"])
Order   = [:constant, :function, :type]
Private = true
```

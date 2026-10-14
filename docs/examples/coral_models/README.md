# ADRIAIndicators with C~scape and CoralPOP output

The data and code for the [Comparing C~scape and CoralPOP](../../src/examples/cscape_coralpop.md)
documentation example, which explains the code step by step.

```
julia --project=. run_all.jl
```

Run this from the folder. It reads the saved output of one default run of each model,
calculates the ADRIAIndicators metrics, prints the values in 2099 and writes two figures to
`output/`.

| Path | Content |
|:-----|:--------|
| `run_all.jl` | The example |
| `data/cscape/` | C~scape output (`.rds`) and the colony diameters and juvenile classes of each type |
| `data/coralpop/` | CoralPOP output (`.npz`) and the run's parameters |
| `data/domain/` | The synthetic one-reef domain: reef areas and the heat stress (DHW) |
| `output/` | The figures used in the documentation |

The data was produced on 2026-10-13 from one default run of each model.

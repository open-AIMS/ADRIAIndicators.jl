# ADRIAIndicators.jl

ADRIAIndicators.jl is a Julia package for summarizing outputs from coral reef ecological
models. It provides a standardized, dependency-free way to transform high-dimensional
model outputs, such as coral abundance by reef, time, species and size class, into
lower-dimensional, interpretable metrics.

Functions range from simple aggregations and unit conversions to indices and
regression-based estimators, covering functional diversity, juvenile abundance, shelter
volume, fish biomass and overall reef condition.

## Installation

ADRIAIndicators.jl requires Julia 1.11 or later.

```julia
using Pkg
Pkg.add("ADRIAIndicators")
```

See [Usage](@ref) to get started.

## The ADRIA ecosystem

ADRIAIndicators.jl began as the metrics module of [ADRIA.jl](https://doi.org/10.5281/zenodo.16900470),
a decision-support platform for evaluating coral reef interventions. It was separated out
so that outputs from different models, or observational data, are summarized with the same
definitions whether or not ADRIA.jl is used.

- [ADRIA.jl](https://doi.org/10.5281/zenodo.16900470): computes its metrics with
  ADRIAIndicators.jl and feeds them into its decision-support algorithms.
- [CoralBlox.jl](https://doi.org/10.5281/zenodo.13251118): the coral model used by ADRIA.jl,
  whose spatial-demographic projections are summarized by ADRIAIndicators.jl.
- [CscapeInterface.jl](https://github.com/open-AIMS/CscapeInterface.jl): uses
  ADRIAIndicators.jl to post-process C-scape simulation results for analysis in ADRIA.jl.
- [ReefGuide](https://doi.org/10.5281/zenodo.20300810): a web-based assessment platform,
  with ADRIAIndicators.jl support planned.

## Available metrics

- [Relative/Absolute/LTMP Cover](@ref "Cover Metrics")
- [Relative/Absolute Shelter Volume](@ref "Metrics")
- [Relative/Absolute Juveniles](@ref "Juvenile Metrics")
- [Juvenile Indicator](@ref "Juvenile Metrics")
- [Coral Diversity](@ref "Metrics")
- [Coral Evenness](@ref "Metrics")
- [Reef Indices](@ref "Reef Indices")

## Citing

If you use ADRIAIndicators.jl, please cite it using the metadata in
[CITATION.cff](https://github.com/open-AIMS/ADRIAIndicators.jl/blob/main/CITATION.cff), or the
Zenodo archive: [10.5281/zenodo.17032475](https://doi.org/10.5281/zenodo.17032475).

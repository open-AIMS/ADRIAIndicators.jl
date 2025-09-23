# ADRIAIndicators.jl

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17032475.svg)](https://doi.org/10.5281/zenodo.17032475)
[![Documentation](https://img.shields.io/badge/docs-dev-blue)](https://open-aims.github.io/ADRIAIndicators.jl/dev/)

ADRIAIndicators provides a set of standard metrics for summarizing the state of reef ecological
model outputs.

## Usage

Each metrics has an option to perform the computation in place.

```julia
using ADRIAIndicators: relative_juveniles, relative_juveniles!

# Raw model coral cover outputs with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations]
raw_model_cover::Array{T,4} = ...
# Juveniles mask with dimensions [sizes]
is_juvenile::Vector{Bool} = ...

# Calculate and allocate new array for metric
rel_juveniles = relative_juveniles(raw_model_cover, is_juvenile)

# Perform the computation inplace.
rel_juveniles = zeros(Float64, n_timesteps, n_locations)
relative_juveniles!(raw_model_cover, is_juvenile, rel_juveniles)
```

## Available Metrics

- Relative/Absolute Cover
- Relative/Absolute Shelter Volume
- Relative/Absolute Juveniles
- Juvenile Indicator
- Coral Diversity
- Coral Evenness

## Reef Indices

The repository also provides implementations of the following indices:

- Reef Condition Index
- Reef Biodiversity Condition Index
- Reef Tourism Index
- Reef Fish Index

They are described in:



## Building Documentation

The documentation is not currently hosted online but can be built as follows.

```bash
cd docs

julia --project=. make.jl
```

finally, open `index.html` in the `docs/build` directory.

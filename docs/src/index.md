# ADRIAIndicators.jl

ADRIAIndicators.jl is a Julia package for summarizing outputs from coral reef ecological
models. Its primary purpose is to provide a standardized and dependency-free package for
transforming high-dimensional model outputs such as coral abundance by reef, time,
species, and size class into lower-dimensional, interpretable metrics.

The package offers a wide range of functions, from simple aggregations and unit
conversions to more complex indices and estimators derived from regression models.
These tools help with the estimation of functional diversity, juvenile abundance,
shelter volume, fish biomass, and overall reef condition, enabling consistent and
comparable analysis across different coral ecology models.


## Usage

```julia
using ADRIAIndicators: relative_juveniles, relative_juveniles!

# Raw model coral cover outputs with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations]
raw_model_cover::Array{T,4} = ...
# Juveniles mask with dimensions [sizes]
is_juvenile::Vector{Bool} = ...
n_timesteps, _, _, n_locations = size(raw_model_cover)

# Calculate and allocate new array for metric
rel_juveniles = relative_juveniles(raw_model_cover, is_juvenile)

# Write the result into provided buffer.
rel_juveniles_out = zeros(Float64, n_timesteps, n_locations)
relative_juveniles!(raw_model_cover, is_juvenile, rel_juveniles_out)
```

Each metric has the option to write the result into an already provided array. This was done
with the intention of writing wrappers in python or R, where memory that is not managed by
the Julia runtime can be written to.

## Available Metrics

- [Relative/Absolute/LTMP Cover](@ref "Cover Metrics")
- [Relative/Absolute Shelter Volume](@ref "Metrics")
- [Relative/Absolute Juveniles](@ref "Juvenile Metrics")
- [Juvenile Indicator](@ref "Juvenile Metrics")
- [Coral Diversity](@ref "Metrics")
- [Coral Evenness](@ref "Metrics")
- [Reef Indices](@ref "Reef Indices")

## Building Documentation

The documentation is hosted online but can also be built locally for offline use or development purposes.

First, ensure you have the project dependencies installed:
```bash
cd docs
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Then, build the documentation:
```bash
julia --project=. make.jl
```

Finally, open `index.html` in the `docs/build` directory to view the documentation.

# ADRIAIndicators.jl

ADRIAIndicators provides a set of standard metrics for summarizing the state of reef ecological
model outputs.

## Usage

Each metric has an option to perform the computation in place for efficiency. The in-place
version of a function is denoted by a `!` at the end of the function name.

```julia
using ADRIAIndicators: relative_juveniles, relative_juveniles!

# Raw model coral cover outputs with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations]
raw_model_cover::Array{T,4} = ...
# Juveniles mask with dimensions [sizes]
is_juvenile::Vector{Bool} = ...
n_timesteps, _, _, n_locations = size(raw_model_cover)

# Calculate and allocate new array for metric
rel_juveniles = relative_juveniles(raw_model_cover, is_juvenile)

# Perform the computation inplace.
rel_juveniles_out = zeros(Float64, n_timesteps, n_locations)
relative_juveniles!(raw_model_cover, is_juvenile, rel_juveniles_out)
```

## Available Metrics

- [Relative/Absolute/LTMP Cover](@ref "Cover Metrics")
- [Relative/Absolute Shelter Volume](@ref "Metrics")
- [Relative/Absolute Juveniles](@ref "Juvenile Metrics")
- [Juvenile Indicator](@ref "Juvenile Metrics")
- [Coral Diversity](@ref "Metrics")
- [Coral Evenness](@ref "Metrics")
- [Reef Indices](@ref "Reef Indices")

## Building Documentation

The documentation is not currently hosted online but can be built locally.

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

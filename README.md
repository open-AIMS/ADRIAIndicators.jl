# ADRIAIndicators.jl

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
relative_juveniles(raw_model_cover, is_juvenile, rel_juveniles)
```

## Available Metrics

- Relative/Absolute Cover
- Relative/Absollute Shelter Volume
- Relative/Absolute Juveniles
- Juvenile Indicator
- Coral Diversity
- Coral Evenness

- Reef Condition Index
- Reef Biodiversity Condition Index
- Reef Tourism Index
- Reef Fish Index

## Building Documentation

The documentation is not currently hosted online but can be built as follows.

```bash
cd docs

julia --project=. make.jl
```

finally, open `index.html` in the `docs/build` directory.

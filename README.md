# ADRIAIndicators.jl

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17032475.svg)](https://doi.org/10.5281/zenodo.17032475)
[![Documentation](https://img.shields.io/badge/docs-dev-blue)](https://open-aims.github.io/ADRIAIndicators.jl/dev/)

ADRIAIndicators provides a set of standard metrics for summarizing the state of reef ecological
model outputs.

## Installation

To install ADRIAIndicators.jl, open a Julia REPL and run:
```julia
using Pkg
Pkg.add("ADRIAIndicators")
```

## Testing

To run tests in a project environment:
```julia
julia> ]test ADRIAIndicators
```

To run tests in a development/testing environment:
```julia
julia> ]test
```

## Usage

Each metric has an option to write the metric into a provided buffer, this version of the function is denoted by a `!` at the end of the function name.

```julia
using ADRIAIndicators

# Create some dummy data
n_timesteps = 10
n_groups = 6
n_sizes = 6
n_locations = 5

# Raw model coral cover outputs with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations]
raw_model_cover = rand(Float64, n_timesteps, n_groups, n_sizes, n_locations);

# Juvenile mask with dimensions [groups ⋅ sizes]. Each row corresponds to a functional group,
# and each column to a size class. Here, the first two size classes are juvenile for all groups.
is_juvenile = repeat([true, true, false, false, false, false]', n_groups, 1);

# Calculate and allocate new array for metric
rel_juveniles = relative_juveniles(raw_model_cover, is_juvenile);
# Output dimensions: [timesteps ⋅ locations] = (10, 5)

# For convenience, if juvenile size classes are the same across all groups,
# a boolean vector can be passed instead:
is_juvenile_vec = [true, true, false, false, false, false];  # [sizes]
rel_juveniles_vec = relative_juveniles(raw_model_cover, is_juvenile_vec);

# Perform the computation and write the metric into a provided buffer.
rel_juveniles_out = zeros(Float64, n_timesteps, n_locations);
relative_juveniles!(raw_model_cover, is_juvenile, rel_juveniles_out);
```

## Available Metrics

- Relative/Absolute/LTMP Cover
- Relative/Absolute Shelter Volume
- Relative/Absolute Juveniles
- Juvenile Indicator
- Coral Diversity
- Coral Evenness
- Reef Condition Index
- Reef Biodiversity Condition Index
- Reef Tourism Index
- Reef Fish Index

## Contributing

Contributions are welcome! Please see the [contributing guidelines](CONTRIBUTING.md) for more information.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Building Documentation

To build documentation locally:

```bash
cd docs

julia --project=. make.jl
```

Then open `index.html` in the `docs/build` directory.

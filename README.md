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
In the ADRIAIndicators.jl testing environment,
```julia
julia> ]test
```

## Usage

Each metric has an option to perform the computation in place. The in-place version of a function is denoted by a `!` at the end of the function name.

```julia
using ADRIAIndicators

# Create some dummy data
n_timesteps = 10
n_groups = 6
n_sizes = 6
n_locations = 5

# Raw model coral cover outputs with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations]
raw_model_cover = rand(Float64, n_timesteps, n_groups, n_sizes, n_locations);

# Juveniles mask with dimensions [sizes]
is_juvenile = [true, true, false, false, false, false];

# Calculate and allocate new array for metric
rel_juveniles = relative_juveniles(raw_model_cover, is_juvenile);

# Perform the computation inplace.
rel_juveniles_out = zeros(Float64, n_timesteps, n_locations);
relative_juveniles!(raw_model_cover, is_juvenile, rel_juveniles_out);
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
*TODO: Add reference to the paper describing the indices.*

## Contributing

Contributions are welcome! Please see the [contributing guidelines](CONTRIBUTING.md) for more information.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use ADRIAIndicators.jl in your research, please cite it as follows:
```bibtex
@misc{Tan_2025,
  author = {Daniel Tan},
  title = {ADRIAIndicators.jl},
  year = {2025},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/open-AIMS/ADRIAIndicators.jl}},
}
```
You can also find the citation information in the `CITATION.cff` file.

## Building Documentation

To build documentation locally:

```bash
cd docs

julia --project=. make.jl
```

Then open `index.html` in the `docs/build` directory.

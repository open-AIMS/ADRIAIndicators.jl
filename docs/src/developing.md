# Developing

## Running tests

From the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Building documentation

The documentation is hosted online, and can also be built locally. The `@example` blocks
are run during the build, so a broken example fails it.

From the repository root, install the dependencies:

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
```

Then build:

```bash
julia --project=docs docs/make.jl
```

The built documentation is `docs/build/index.html`.

## Cookbook figures

The cookbook pages use plain code blocks (they are not run by the documentation build), so
their figures are committed to `docs/src/assets/cookbook/`.

- The CoralBlox.jl figures are regenerated from the page's own code with
  `julia --project=docs/figures docs/figures/coralblox_figures.jl` (the environment is
  instantiated with `julia --project=docs/figures -e 'using Pkg; Pkg.instantiate()'`).
- The C~scape and CoralPop figures are written by `examples/coral_models/run_all.jl` to its
  `output/` folder, from where they are copied.


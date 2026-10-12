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

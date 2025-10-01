# Contributing to ADRIAIndicators.jl

Thank you for considering contributing to ADRIAIndicators.jl.

## Reporting Bugs

Bugs are tracked as [GitHub issues](https://github.com/open-AIMS/ADRIAIndicators.jl/issues). Please include the following information in your bug report:

*   A clear and descriptive title.
*   Steps to reproduce the behavior.
*   A minimal, reproducible example.
*   An explanation of what you expected to happen and what actually happened.
*   Your operating system, Julia version, and ADRIAIndicators.jl version.

## Pull Request Process

1.  Ensure that all tests pass.
2.  Format your code with `JuliaFormatter.jl` before submitting.

## Development Setup

1.  Fork the repository and clone it to your local machine.
2.  Install the development dependencies:
    ```julia
    using Pkg
    Pkg.activate(".")
    Pkg.instantiate()
    ```
3.  Run the tests:
    ```julia
    using Pkg
    Pkg.test("ADRIAIndicators")
    ```

## Code Style

This project follows the [Blue Style Guide](https://github.com/invenia/BlueStyle) for Julia code. We use [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl) to enforce this style. Before submitting a pull request, please format your code by running:
```julia
using JuliaFormatter
format(".")
```

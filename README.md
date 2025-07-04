# ReefMetrics.jl

## Build Instructions

```bash
cd path/to/ReefMetrics.jl

cd build

julia --project=. build.jl
```

## Supported Metrics

**Reef State Metrics**

- Relative Cover
- Total Absolute Cover
- Relative Taxa Cover
- Relative Loc Taxa Cover
- Relative Juveniles
- Absolute Juveniles
- Coral Diversity
- Coral Evenness
- Juveniles Indicator
- Absolute Shelter Volume
- Relative Shelter Volume

**Biodiversity Metrics**

- Reef Biodiversity Condition Index

**Economic Metrics**

- Reef Tourism Index
- Reef Condition Index
- Reef Fish Index

**Others**

- Simpsons Diversity

## Function inputs

Function accept inputs in the following order,
1. Size of dimensions, with order prioritised by
    1. Time
    2. Groups
    3. Size
    4. Locations
2. Input Arrays with documented dimensions,
3. Non Array inputs
4. Output array buffer with documented dimensions

```julia
Base.@ccallable function metric(
    n_timesteps::Cint,
    n_groups::Cint,
    n_size::Cint,
    n_locations::Cint, # size of dimensions
    input_array1::Ptr{Float64},
    ...
    non_array_input1::Cdouble,
    non_array_input2::Clong,
    preallocated_output::Ptr{Float64}
)::Cvoid
    # Wrap input matrices
    inputs_array1::Array{Float64, n} = unsafe_wrap(Array, input_array1, (dims...))
    ...
    # Call julia implementations
    metric(input_array1 ..., non_array_inputs..., preallocated_output...)

    return nothing
end
```

### Examples

```julia
"""
    reef_condition_index(n_timesteps::Cint, n_locs::Cint, relative_cover::Ptr{Float64}, evenness::Ptr{Float64}, relative_shelter_volume::Ptr{Float64}, juvenile_indicator::Ptr{Float64}, threshold::Cint, output_rci::Ptr{Float64})

# Arguments:
- `n_timesteps`: Number of timesteps in input data arrays.
- `n_locs`: Number of locations
- `relative_cover`: Relative cover with dimensions [timesteps ⋅ locations]
- `evenness`: Evenness with dimensions [timesteps ⋅ locations]
- `relative_shelter_volume`: RSV with dimensions [timesteps ⋅ locations]
- `juvenile_indicator`: Juveniles indicator with dimensions [timesteps ⋅ locations]
- `threshold`: Conditions thresholds
- `output_rci`: Output array buffer of dimensions [timesteps ⋅ locations]
"""
Base.@ccallable function reef_condition_index(
    n_timesteps::Cint,
    n_locs::Cint,
    relative_cover::Ptr{Float64},
    evenness::Ptr{Float64},
    relative_shelter_volume::Ptr{Float64},
    juvenile_indicator::Ptr{Float64},
    threshold::Cint,
    output_rci::Ptr{Float64}
)::Cvoid
    return nothing
end
```

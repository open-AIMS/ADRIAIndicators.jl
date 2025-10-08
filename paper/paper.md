---
title: 'ADRIAIndicators: a Julia package for summarizing reef ecological model outputs'
tags:
  - Julia
  - ecology
  - coral reefs
  - metrics
authors:
  - name: Daniel Tan
    orcid: 0009-0004-8696-0631
    affiliation: 1
  - name: Takuya Iwanaga
    orcid: 0000-0001-8173-0870
    affiliation: 1
affiliations:
 - name: Australian Institute of Marine Science
   index: 1
date: 29 September 2025
bibliography: paper.bib
---

# Summary

ADRIAIndicators.jl is a Julia package for analyzing outputs from coral reef ecological
models. Its primary purpose is to provide a standardized and dependency-free toolkit for
transforming high-dimensional model outputs such as coral abundance by reef, time,
species, and size class into lower-dimensional, interpretable metrics. The package offers a
wide range of functions, from simple aggregations and unit conversions to more complex
indices and estimators derived from regression models. These tools help the
assessment of species diversity, juvenile abundance, shelter volume, fish biomass, and
overall reef condition, enabling consistent and comparable analysis across different coral
ecology models such as CoralBlox [@CoralBlox], C-Scape [@CScape], and ReefMod [@ReefMod].

# Statement of Need

Ecological models of coral reefs produce large amounts of high-dimensional data. There is a
need for standardized tools to summarize and analyze these model outputs to facilitate
comparison between different models and scenarios. ADRIAIndicators.jl provides a set of
standard metrics that can be used to summarize the state of reef ecological model outputs
that previously existed within the ADRIA.jl package but was being reproduced for other coral
ecology models. The package is written in Julia, a high-level, high-performance programming language for
technical computing. The package is designed to be easy to use, and
provides an in-place option for all metrics for any eventual wrappers that may be
implemented in python and R.

## Available Metrics

The functions implemented in ADRIAIndicators.jl are classified into three categories:
Metrics, Aggregations, and Conversions. Aggregations are convenience functions for
reducing the dimensionality of data by summarizing arrays. Conversions handle
transformations between different units or representations of coral cover.
Metrics are functions that derive higher-level, interpretable indicators from the raw model
data, such as coral diversity, shelter volume, and composite indices for reef health.

| **Metric Name**                                  | **Type**          | **Reference**   |
|--------------------------------------------------|-------------------|-----------------|
| Absolute Shelter Volume                          | Metric            | [@URBINABARRETO2021107151; @ASTON_STRUCTURAL]|
| Relative Shelter Volume                          | Metric            | -               |
| Coral Diversity                                  | Metric            | [@CoralDiversity]|
| Coral Evenness                                   | Metric            | -               |
| Reef Condition Index                             | Metric            |[@ReefConditionIndex]|
| Reef Tourism Index                               | Metric            |               |
| Reef Biodiversity Condition Index                | Metric            |               |
| Reef Fish Index                                  | Metric            | [@ReefFishIndex]               |
| Relative Cover                                   | Aggregation       |                 |
| Relative Location Taxonomy Cover                 | Aggregation       |                 |
| Relative Taxonomy Cover                          | Aggregation       |                 |
| LTMP Cover                                       | Aggregation       |                 |
| LTMP Location Taxonomy Cover                     | Aggregation       |                 |
| LTMP Taxonomy Cover                              | Aggregation       |                 |
| Relative Juveniles                               | Aggregation       |                 |
| Relative Location Taxonomy Juveniles             | Aggregation       |                 |
| Relative Taxonomy Juveniles                      | Aggregation       |                 |
| Relative Habitable Cover to Reef Cover           | Conversion        |                 |
| Reef Cover to Relative Habitable Cover           | Conversion        |                 |
                                                   

## Usage

The order of dimensions is always the same in ADRIAIndicators.jl,

1. Time
2. Groups
3. Sizes
4. Locations
5. Scenarios

If a dimension is missing then the order remains the same however the missing dimensions is
excluded. Furthermore, all metrics have an option to provide a buffer as input in the cases
where one wants to write the metric into an existing array or sub-array. This implementation
is relied upon by functions that allocate the returned array as-well and was chosen to
account for any any decisions in the future where another language may wrap this library and
need to pass memory that is not managed by Julia.

```julia
using ADRIAIndicators

# Create some dummy model output data
n_timesteps = 75
n_groups = 5
n_sizes = 7
n_locations = 3806

# Raw model coral cover outputs with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations]
raw_model_cover = rand(Float64, n_timesteps, n_groups, n_sizes, n_locations);

# Juveniles mask with dimensions [sizes]
is_juvenile = [true, true, false, false, false, false, false];

# Calculate and allocate new array for metric
rel_juveniles = relative_juveniles(raw_model_cover, is_juvenile);

# Users can provide output buffers if it is more convenient for them
rel_juveniles_out = zeros(Float64, n_timesteps, n_locations);
relative_juveniles!(raw_model_cover, is_juvenile, rel_juveniles_out);
```

# Acknowledgements

This package originally existed as a module of the ADRIA.jl platform and is managed by the
Australian Institute of Marine Science (AIMS) Decision Support team.

# References

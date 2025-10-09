---
title: 'ADRIAIndicators.jl: a Julia package for summarizing reef ecological model outputs'
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
  - name: Rose Crocker
    orcid: 0009-0007-3586-705X
    affiliation: 1
  - name: Ken Anthony
    orcid:
    affiliation: 3
  - name: Pedro Ribeiro de Almeida
    orcid: 0009-0007-4814-3261
    affiliation: 1
  - name: Arne Adam
    orcid: 0000-0002-2960-7880
    affiliation: 2
  - name: Ryan F. Heneghan
    orcid: 0000-0001-7626-1248
    affiliation: 4
  - name: Michael McWilliam
    orcid: 0000-0001-5748-0859
    affiliation: 5
  - name: Morgan Pratchett
    orcid: 0000-0002-1862-8459
    affiliation: 5
  - name: Vanessa Haller-Bull
    orcid: 0000-0002-3919-7053
    affiliation: 1
  - name: Yves-Marie Bozec
    orcid: 0000-0002-7190-5187
    affiliation: 2
  - name: Anna K. Cresswell
    orcid: 0000-0001-6740-9052
    affiliation: 1
  - name: Juan Carlos Ortiz
    orcid:
    affiliation: 1
affiliations:
 - name: Australian Institute of Marine Science, Townsville, Queensland, Australia
   index: 1
 - name: School of the Environment, University of Queensland, St Lucia, Queensland, Australia
   index: 2
 - name: Nature Assets Consulting
   index: 3
 - name: School of Environment and Science, Griffith University, Nathan 4111, QLD Australia
   index: 4
 - name: College of Science and Engineering, James Cook University, Townsville, QLD 4811, Australia
   index: 5
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
ecology models such as CoralBlox [@CoralBlox], C-Scape [@CScape], ReefMod [@ReefMod], and
CoCoNet *[citation]*.

# Statement of Need

Models of coral reef ecosystems can produce large volumes of high-dimensional data. There
is a need for standardized tools to summarize and analyze these model outputs to facilitate
inter-model comparison of environmental projections and implications of intervention
activities. ADRIAIndicators.jl provides a set of standard indicator metrics that can be used
to summarize the state of reef ecological model outputs that previously existed within the
ADRIA.jl Decision Support package *[citation]* but was being reproduced for other coral
ecology models. ADRIAIndicators.jl is written in Julia *[citation to https://doi.org/10.1137/141000671]*,
a high-level, high-performance programming language for technical computing. This package
is designed to be easy to use, and provides an in-place option for all metrics for any
eventual wrappers that may be implemented in other languages such as Python and R.

## Available Metrics

The indicators implemented in ADRIAIndicators.jl are classified into three categories:
Metrics, Aggregations, and Conversions. Aggregations are convenience methods for
reducing the dimensionality of data by summarizing arrays. Conversions handle
transformations between different units or representations of coral cover.
Metrics derive higher-level, interpretable indicators from the raw model
data, such as coral diversity, shelter volume, and composite indices for reef health.

| **Metric Name**                                  | **Type**          | **Reference**   |
|--------------------------------------------------|-------------------|-----------------|
| Absolute Shelter Volume                          | Metric            | [@URBINABARRETO2021107151; @ASTON_STRUCTURAL]|
| Relative Shelter Volume                          | Metric            | -               |
| Coral Diversity                                  | Metric            | [@CoralDiversity]|
| Coral Evenness                                   | Metric            | -               |
| Reef Condition Index                             | Metric            |[@ReefConditionIndex]|
| Reef Tourism Index                               | Metric            |               |
| Reef Biodiversity Condition Index                | Metric            | *[citation MW paper]* |
| Reef Fish Index                                  | Metric            | [@ReefFishIndex] |
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

AIMS acknowledges the Traditional Owners of the Land and Sea Country on which these
indicators were developed. This package originally existed as a module of the ADRIA.jl
platform and is maintained by the Australian Institute of Marine Science (AIMS) Decision
Support team. Work was conducted and funded as part of the Reef Restoration and Adaptation
Program (RRAP) and includes additional work funded by the Great Barrier Reef Foundation
(GBRF).

# References

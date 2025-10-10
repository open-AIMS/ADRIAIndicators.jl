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
  - name: Pedro Ribeiro de Almeida
    orcid: 0009-0007-4814-3261
    affiliation: 1
  - name: Ken Anthony
    orcid:
    affiliation: 2
  - name: Arne A. S. Adam
    orcid: 0000-0002-2960-7880
    affiliation: 3
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
    affiliation: 3
  - name: Anna K. Cresswell
    orcid: 0000-0001-6740-9052
    affiliation: 1
  - name: Juan Carlos Ortiz
    orcid:
    affiliation: 1
  - name: Scott Condie
    orcid: 0000-0002-5943-014X
    affiliation: 6
affiliations:
 - name: Australian Institute of Marine Science, Townsville, Queensland, Australia
   index: 1
 - name: Nature Assets Consulting
   index: 2
 - name: School of the Environment, University of Queensland, St Lucia, Queensland, Australia
   index: 3
 - name: School of Environment and Science, Griffith University, Nathan 4111, QLD Australia
   index: 4
 - name: College of Science and Engineering, James Cook University, Townsville, QLD 4811, Australia
   index: 5
 - name: CSIRO Environment, Hobart, TAS, Australia
   index: 6
date: 29 September 2025
bibliography: paper.bib
---

# Summary

ADRIAIndicators.jl is a Julia package for analyzing outputs from coral reef ecological
models. Its primary purpose is to provide a standardized and dependency-free toolkit for
transforming high-dimensional model outputs such as coral abundance by reef, time,
species, and size class into lower-dimensional, interpretable metrics. The package offers a
wide range of functions, from simple aggregations and unit conversions to more complex
indices and estimators derived from regression models. These tools help with the
estimation of functional diversity, juvenile abundance, shelter volume, fish biomass, and
overall reef condition, enabling consistent and comparable analysis across different coral
ecology models such as CoralBlox [@CoralBlox], C~Scape [@CScape], ReefMod [@ReefMod], and
CoCoNet [@CoCoNet].

# Statement of Need

Models of coral reef ecosystems often produce large volumes of high-dimensional data. There
is a need for standardized tools to summarize and analyze these model outputs to facilitate
inter-model comparison of environmental projections and communicate results to managers and
stakeholders. ADRIAIndicators.jl provides a set of standard indicator metrics that can be
used to summarize reef state in ecological model outputs that previously existed within
the ADRIA.jl Decision Support package [@ADRIA] but were being reproduced in many workflows
that did not use ADRIA.jl.

ADRIAIndicators.jl is written in Julia [@Julia], a high-level, high-performance programming
language for technical computing. This package is designed to be easy to use, and provides
an in-place option for all metrics for any eventual wrappers that may be implemented in
other languages such as Python and R. Such wrappers could be developed leveraging Julia's
support for language interoperability and compilation capabilities (such as those provided
by JuliaC.jl).

## Available Indicators

The indicators implemented in ADRIAIndicators.jl are classified into three categories:
Aggregations, Conversions, and Metrics. Aggregations are convenience methods for
reducing the dimensionality of data by summarizing arrays. Conversions handle
transformations between different units or representations of coral cover.
Metrics derive higher-level, interpretable indicators from the raw model
data, such as coral diversity, shelter volume, and composite indices for reef health.

| **Metric Name**                                  | **Type**          | **Reference**   |
|--------------------------------------------------|-------------------|-----------------|
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
| Absolute Shelter Volume                          | Metric            | [@URBINABARRETO2021107151; @ASTON_STRUCTURAL]|
| Relative Shelter Volume                          | Metric            | -               |
| Coral Diversity                                  | Metric            | [@CoralDiversity]|
| Coral Evenness                                   | Metric            | -               |
| Reef Condition Index                             | Metric            |[@ReefConditionIndex]|
| Reef Tourism Index                               | Metric            |               |
| Reef Biodiversity Condition Index                | Metric            | *[citation MW paper]* |
| Reef Fish Index                                  | Metric            | [@ReefFishIndex] |

*A dash (-) in the 'Reference' column indicates the reference is the same as the entry directly above it.*

### Indicator Summaries

Each indicator is briefly summarized below. Full implementation details are found in the
documentation, including descriptions of their mathematical formulations where appropriate.

#### Coral Cover

Estimates of coral cover are provided in both **Relative** and **Absolute** forms and
estimated for each location by summing over functional groups and their size classes. As
the indicators are agnostic to the spatial scale being assessed, the term "location" is used
to convey an arbitrary unit of analysis. For example, a location could be a representative
reef, site within a reef, a transect, or patch/plot of reef.

Relative cover is calculated to be *relative* to the location's coral habitable area. The
meaning of "habitable" area is subject to much debate, however it can be construed as being
representative of the area of hard substrate or a location's carrying capacity. The term
*LTMP cover* is used to convey that cover estimates are made relative to estimates of the
total reef area, inclusive of reef areas where corals are unable to settle. As such,
LTMP cover estimates may never reach 100%. This approach is more in line with the values
reported by the Long-Term Monitoring Program (LTMP). The *Absolute* form provides estimates
of the area of coral cover expressed in SI units (typically m²).

Relevant indicators:

-   **Relative Cover**: Calculates the total relative cover per location by summing over functional groups and size classes.
-   **Relative Location Taxonomy Cover**: Calculates the relative cover for each location and functional group by aggregating size classes.
-   **Relative Taxonomy Cover**: Provides an indication of the coral cover decomposed by functional group by aggregating their size classes for all locations.
-   **LTMP Cover**: Calculates the coral cover for each location relative to estimated total reef area. More comparable to the values reported by the Long Term Monitoring Program (LTMP).
-   **LTMP Location Taxonomy Cover**: As above, but decomposes the coral cover estimates to each functional group by location.
-   **LTMP Taxonomy Cover**: As above, but providing total values per functional group.
-   **Relative Juveniles**: Calculates the relative coral cover provided by juvenile corals. User indicates which size classes are construed to be "juvenile".
-   **Relative Location Taxonomy Juveniles**: As above, but for each location and functional group.
-   **Relative Taxonomy Juveniles**: As above, but summed across all locations.

#### Shelter Volume

Calculates the volume of shelter provided by the given coral cover. In typical use, values
are indicative of the modelled *live* coral population, however it is noted that non-living
substrate may also provide some form of shelter.

Relevant indicators:

-   **Absolute Shelter Volume**: Calculates the absolute shelter volume (in m³) provided by corals.
-   **Relative Shelter Volume**: Calculates the relative shelter volume, expressed as a proportion of the theoretical maximum shelter volume for a given area.

#### Diversity and Evenness

These indicators provide estimates of the diversity and evenness of coral functional groups.

Relevant indicators:

-   **Coral Diversity**: Calculates coral diversity at each location using the Simpson's Diversity Index, which accounts for the number and relative abundance of coral functional groups.
-   **Coral Evenness**: Calculates the evenness of coral functional groups at each location using the Inverse Simpson's Index, indicating how similar in abundance the different functional groups are.

#### Condition Indices

Composite indices (or meta-metrics; metrics of metrics) that provide a single value
indication of reef condition(s).

-   **Reef Condition Index**: A categorical index (from 'Very Poor' to 'Very Good') that assesses overall reef health based on coral cover, shelter volume, juvenile abundance, and rubble cover.
-   **Reef Tourism Index**: A continuous index fitted with a linear regression model to assess reef health for tourism purposes, based on metrics like coral cover, evenness, shelter volume, and juvenile abundance.
-   **Reef Fish Index**: An index that estimates fish biomass based on a relationship between coral cover and structural complexity.
-   **Reef Biodiversity Condition Index**: An index calculated as the average of relative coral cover, coral diversity, and relative shelter volume to represent the state of reef biodiversity.

#### Conversions

These are convenience/helper methods to aid in data transformations to various forms.

Relevant methods:

-   **Relative Habitable Cover to Reef Cover**: Converts relative coral cover (proportion of habitable area) to LTMP cover (proportion of total reef area).
-   **Reef Cover to Relative Habitable Cover**: Converts LTMP cover to relative coral cover.


## Usage

The order of dimensions is always the same in ADRIAIndicators.jl,

1. Time
2. Groups
3. Sizes
4. Locations
5. Scenarios

If a dimension is missing then the order remains the same however the missing dimensions are
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

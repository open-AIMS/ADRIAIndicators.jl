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
affiliations:
 - name: Australian Institute of Marine Science
   index: 1
date: 29 September 2025
bibliography: paper.bib
---

# Summary

ADRIAIndicators is a Julia package designed to transform high-dimensional data from coral
ecological models such as CoralBlox, C-Scape, and ReefMod into meaningful metrics [@ReefMod; @CScape].
It provides a suite of functions to process modelled coral abundance by reef, time step,
species, and size class and transforms them into lower-dimensional, interpretable
indicators. The package includes metrics for assessing species diversity, juvenile abundance,
and shelter volume, alongside tools for data aggregation.

# Statement of Need

Ecological models of coral reefs are complex and produce large amounts of data. There is a
need for standardized tools to summarize and analyze these model outputs to facilitate
comparison between different models and scenarios. ADRIAIndicators provides a set of
standard metrics that can be used to summarize the state of reef ecological model outputs.
The package is written in Julia, a high-level, high-performance programming language for
technical computing. The package is designed to be easy to use, and it
provides an in-place option for all metrics for any eventual wrappers that may be
implemented in python and R.

# Available Metrics

Implemented functions are classified by metrics, aggregations and conversions. Aggregations
are simply convenience functions used to reduce the dimensionality of the input array and
summarise high dimensional arrays.

The conversion of relative habitable cover to reef cover
is important because of two differing types of coral cover that are typically reported
between coral ecologists. Where relative habitable cover refers to the proportion of area
occupied by coral relative to area that is physically able to be occupied by corals,
where as reef cover refers to the proportion of `reef` that is occupied by corals regardless
of whether some area is not able to lived on by coral.

The reef condition index, reef fish index, and reef tourism index are indices based on
expert solicitation of what constitutes a healthy reef and regressions [@ReefFishIndex].


| **Type**        | **Metric Name**                                  | **Reference**   |
|-----------------|--------------------------------------------------|-----------------|
| Metric            | Absolute Shelter Volume                          | [@URBINABARRETO2021107151; @ASTON_STRUCTURAL]|
| Metric            | Relative Shelter Volume                          | -               |
| Metric            | Coral Diversity                                  | [@CoralDiversity]|
| Metric            | Coral Evenness                                   | -               |
| Metric            | Reef Condition Index                             | [@ReefFishIndex]|
| Metric            | Reduced Reef Condition Index                     | -               |
| Metric            | Reef Tourism Index                               | -               |
| Metric            | Reduced Reef Tourism Index                       | -               |
| Metric            | Reef Fish Index                                  | -               |
| Metric            | Reef Biodiversity Condition Index                | -               |
| Aggregation       | Relative Cover                                   |                 |
| Aggregation       | Relative Location Taxonomy Cover                 |                 |
| Aggregation       | Relative Taxonomy Cover                          |                 |
| Aggregation       | Relative Juveniles                               |                 |
| Aggregation       | Relative Location Taxonomy Juveniles             |                 |
| Aggregation       | Relative Taxonomy Juveniles                      |                 |
| Conversion       | Relative Habitable Cover to Reef Cover           |                 |
| Conversion       | Reef Cover to Relative Habitable Cover           |                 |


# Acknowledgements

# References

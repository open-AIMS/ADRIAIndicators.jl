#=
pipeline.jl -- a model's coral cover  ->  ADRIAIndicators  ->  indicators per location and per domain.

Nothing in here is specific to C~scape or CoralPop: any model whose output can be put into the "standard
layout" (see `check_model`) goes through exactly the same code. `adapters.jl` does that for the two example
models; `run_all.jl` runs everything and draws the figures; README.md explains every step.
=#

import ADRIAIndicators as AI
using Statistics

# ===== constants both models share =====

## Reference colony that relative shelter volume is scaled against: one tabular Acropora of 95 cm
## diameter per m2, given as (diameter in cm, intercept b, coefficient a) of log(S) = b + a log(x).
## This is ADRIA's convention; it is what makes shelter volume comparable between models.
const SV_REFERENCE = (95.0, -8.95, 2.80)

## Juvenile density (colonies per m2) at which the juvenile indicator reaches 1 (ADRIA's value).
const MAX_JUVENILE_DENSITY = 15.0

# ===== check that a model is in the standard layout =====

"""
    check_model(m)

The standard layout, one NamedTuple per model (T timesteps, G functional groups, S size classes,
L locations, N runs; N = 1 in this example):

| field         | shape           | meaning                                                                |
|---------------|-----------------|------------------------------------------------------------------------|
| `cover`       | `[T, G, S, L, N]` | cover of each group / size class as a FRACTION of the location's habitable (k) area |
| `k_area`      | `[L]`           | habitable area of each location (m2)                                   |
| `reef_area`   | `[L]`           | total area of each location (m2)                                       |
| `diam_cm`     | `[G, S, N]`     | mean colony diameter of each group / size class (cm)                   |
| `pap`         | `[G, S, 2]`     | planar-area parameters (intercept, coefficient)                        |
| `is_juvenile` | `[G, S]`        | which group / size classes count as juvenile                           |

plus `label`, `years`, `groups` and `locations`. Stops with a readable message if a shape or range is wrong.
"""
function check_model(m)
    T, G, S, L, N = size(m.cover)
    length(m.years) == T || error("$(m.label): $(length(m.years)) years but cover has $T timesteps")
    length(m.groups) == G || error("$(m.label): $(length(m.groups)) group names but cover has $G groups")
    size(m.diam_cm) == (G, S, N) || error("$(m.label): diam_cm is $(size(m.diam_cm)), expected $((G, S, N))")
    size(m.pap) == (G, S, 2) || error("$(m.label): pap is $(size(m.pap)), expected $((G, S, 2))")
    size(m.is_juvenile) == (G, S) || error("$(m.label): is_juvenile is $(size(m.is_juvenile)), expected $((G, S))")
    length(m.k_area) == L && length(m.reef_area) == L || error("$(m.label): k_area / reef_area must have $L entries")
    all(0 .<= m.cover .<= 1) || error("$(m.label): cover must be a fraction of habitable area in [0, 1]")
    all(m.k_area .<= m.reef_area) || error("$(m.label): habitable area exceeds reef area")
    return nothing
end

# ===== the indicators, per location and run =====

## (key, title shown on the figure)
const METRICS = [
    (:relative_cover, "Relative cover"),
    (:relative_juveniles, "Relative juveniles"),
    (:juvenile_indicator, "Juvenile indicator"),
    (:relative_shelter_volume, "Relative shelter volume"),
    (:absolute_shelter_volume, "Absolute shelter volume (m³)"),
    (:coral_evenness, "Coral evenness"),
    (:coral_diversity, "Coral diversity"),
    (:reef_condition_index, "Reef Condition Index"),
    (:reef_tourism_index, "Reef Tourism Index"),
    (:reef_fish_index, "Reef Fish Index"),
    (:reef_biodiversity_condition_index, "Reef Biodiversity Condition Index"),
]

"""
    compute_indicators(m) -> Dict{Symbol,Array{Float64,3}}

Every indicator of `METRICS` as a `[timesteps, locations, runs]` array, using only ADRIAIndicators.
ADRIAIndicators works on one run at a time (`[T, G, S, L]`), so this loops over the runs.
"""
function compute_indicators(m)
    T, G, S, L, N = size(m.cover)
    out = Dict{Symbol,Array{Float64,3}}(k => zeros(T, L, N) for (k, _) in METRICS)

    for n in 1:N
        cover = m.cover[:, :, :, :, n]                                  # [T, G, S, L]
        diam_cm = m.diam_cm[:, :, n]                                    # [G, S]

        ## -- cover --
        rc = AI.relative_cover(cover)                                   # [T, L]  all groups and sizes
        taxa = AI.relative_loc_taxa_cover(cover)                        # [T, G, L]  per group
        ltmp = AI.ltmp_cover(cover, m.k_area, m.reef_area)              # [T, L]  relative to REEF area

        ## -- shelter volume: ADRIAIndicators returns it per group / size class; add them up and,
        ##    for the relative version, clamp to [0, 1] (the metric's definition) --
        rsv = AI.relative_shelter_volume(cover, diam_cm, m.pap, m.k_area, SV_REFERENCE)
        rsv = clamp.(dropdims(sum(rsv; dims=(2, 3)); dims=(2, 3)), 0.0, 1.0)             # [T, L]
        asv = AI.absolute_shelter_volume(cover, diam_cm, m.pap, m.k_area)
        asv = dropdims(sum(asv; dims=(2, 3)); dims=(2, 3))                               # [T, L], m3

        ## -- juveniles. NOTE the unit: the indicator needs diameters in METRES, shelter volume in CM --
        rjuv = AI.relative_juveniles(cover, m.is_juvenile)              # [T, L]
        jind = AI.juvenile_indicator(cover, m.is_juvenile, m.k_area, diam_cm ./ 100.0, MAX_JUVENILE_DENSITY)

        ## -- diversity --
        even = AI.coral_evenness(taxa)                                  # [T, L]
        div = AI.coral_diversity(taxa)                                  # [T, L]

        ## -- indices built from the metrics above (no rubble: neither model has it) --
        rci = AI.reef_condition_index(ltmp, rsv, jind)
        rti = AI.reef_tourism_index_no_rubble(ltmp, even, rsv, rjuv)
        rfi = AI.reef_fish_index(rc)
        rbci = AI.reef_biodiversity_condition_index(rc, div, rsv)

        for (k, v) in (:relative_cover => rc, :relative_juveniles => rjuv, :juvenile_indicator => jind,
            :relative_shelter_volume => rsv, :absolute_shelter_volume => asv, :coral_evenness => even,
            :coral_diversity => div, :reef_condition_index => rci, :reef_tourism_index => rti,
            :reef_fish_index => rfi, :reef_biodiversity_condition_index => rbci)
            out[k][:, :, n] .= v
        end
    end
    return out
end

# ===== from locations to the domain =====

"""
    domain_metrics(m, ind) -> Dict{Symbol,Matrix{Float64}}

Collapse the location dimension of every indicator: the k-area weighted mean (each location counts in
proportion to its habitable area) -- or, for the one absolute quantity (shelter volume in m3), the sum.
Result: `[timesteps, runs]`. With a single location both are the location's own value.
"""
function domain_metrics(m, ind)
    return Dict(k => k === :absolute_shelter_volume ?
                     AI.scenario_metric(v, m.k_area, 2; is_relative=false, return_relative=false) :
                     AI.scenario_metric(v, m.k_area, 2)
                for (k, v) in ind)
end

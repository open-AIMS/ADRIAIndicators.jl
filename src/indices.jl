"""
    reef_condition_index!(ltmp_cover::AbstractArray{T,2}, relative_shelter_volume::AbstractArray{T,2}, juvenile_indicator::AbstractArray{T,2}, out_rci::AbstractArray{T,2})::Nothing where {T<:AbstractFloat}

Calculate the Reef Condition Index (RCI).

This method uses three inputs: LTMP coral cover, relative shelter volume, and juvenile
indicator. This is implementation is limited to using outputs that can be provided by
ecological models that only model corals.

# Arguments
- `ltmp_cover` : LTMP Coral Cover with dimensions [timesteps ⋅ locations].
- `relative_shelter_volume` : Relative shelter volume with dimensions [timesteps ⋅ locations].
- `juvenile_indicator`: Juvenile Indicator with dimensions [timesteps ⋅ locations].
- `out_rci` : Output RCI buffer with dimensions [timesteps ⋅ locations].

# References
1. Ryan F. Heneghan, Gabriela Scheufele, Yves-Marie Bozec et al. A framework to inform
    economic valuation of non-use benefits from coral-reef intervention efforts, 02 October
    2025, PREPRINT (Version 1) available at Research Square
    [https://doi.org/10.21203/rs.3.rs-7644150/v1]
"""
function reef_condition_index!(
    ltmp_cover::AbstractArray{<:AbstractFloat,2},
    relative_shelter_volume::AbstractArray{<:AbstractFloat,2},
    juvenile_indicator::AbstractArray{<:AbstractFloat,2},
    out_rci::AbstractArray{<:AbstractFloat,2};
)::Nothing
    RCI_CRIT = eltype(ltmp_cover)[
        0.05 0.15 0.25 0.35 0.45;  # relative cover thresholds
        0.15 0.25 0.30 0.35 0.45;  # shelter volume thresholds
        0.15 0.20 0.25 0.30 0.35   # coral juveniles thresholds
    ]

    CRIT_VAL = eltype(ltmp_cover)[0.1, 0.3, 0.5, 0.7, 0.9]

    CRITERIA_THRESHOLD = 2

    out_rci .= 0.1
    metrics_met = zeros(eltype(ltmp_cover), size(ltmp_cover)...)
    pass_mask = falses(size(ltmp_cover)...)
    for i in 1:size(CRIT_VAL, 1)
        @. metrics_met = (ltmp_cover >= RCI_CRIT[1, i]) +
                         (relative_shelter_volume >= RCI_CRIT[2, i]) +
                         (juvenile_indicator >= RCI_CRIT[3, i])
        @. pass_mask = metrics_met >= CRITERIA_THRESHOLD

        out_rci[pass_mask] .= CRIT_VAL[i]
    end


    return nothing
end

"""
    reef_condition_index(ltmp_cover::AbstractArray{T,2}, relative_shelter_volume::AbstractArray{T,2}, juvenile_indicator::AbstractArray{T,2})::AbstractArray{T,2} where {T<:AbstractFloat}

Calculate the Reef Condition Index (RCI).

This method uses three inputs: LTMP cover, relative shelter volume, and juvenile indicator.

The RCI is a categorical index assessing the overall health and condition of a reef location
based on three key ecological metrics. The index assigns a discrete score (0.1, 0.3, 0.5, 0.7, or 0.9)
representing categories from "Very Poor" to "Very Good".

For each input there are five levels of condition ranging from very poor to very good. Then
the location is assigned a condition of very poor to very good if that location meets at
least two of the metrics condition criteria. The condition level is then assigned a
numerical value based on its categorisation.

# Arguments
- `ltmp_cover` : Relative Coral Cover with dimensions [timesteps ⋅ locations].
- `relative_shelter_volume` : Relative shelter volume with dimensions [timesteps ⋅ locations].
- `juvenile_indicator`: Relative juvenile cover with dimensions [timesteps ⋅ locations].

# Returns
Output RCI buffer with dimensions [timesteps ⋅ locations].

# References
1. Ryan F. Heneghan, Gabriela Scheufele, Yves-Marie Bozec et al. A framework to inform
    economic valuation of non-use benefits from coral-reef intervention efforts, 02 October
    2025, PREPRINT (Version 1) available at Research Square
    [https://doi.org/10.21203/rs.3.rs-7644150/v1]
"""
function reef_condition_index(
    ltmp_cover::AbstractArray{<:AbstractFloat,2},
    relative_shelter_volume::AbstractArray{<:AbstractFloat,2},
    juvenile_indicator::AbstractArray{<:AbstractFloat,2}
)::AbstractArray{<:AbstractFloat,2}
    lc_size = size(ltmp_cover)
    if (lc_size != size(relative_shelter_volume)) || (lc_size != size(juvenile_indicator))
        throw(DimensionMismatch("All input metric arrays must have the same dimensions."))
    end

    out_rci = zeros(eltype(ltmp_cover), size(ltmp_cover))
    reef_condition_index!(
        ltmp_cover,
        relative_shelter_volume,
        juvenile_indicator,
        out_rci
    )

    return out_rci
end


"""
    reef_condition_index!(ltmp_cover::AbstractArray{T,2}, relative_shelter_volume::AbstractArray{T,2}, juvenile_indicator::AbstractArray{T,2}, rubble::AbstractArray{T,2}, out_rci::AbstractArray{T,2})::Nothing where {T<:AbstractFloat}

Calculate the Reef Condition Index (RCI).

This method uses four inputs: LTMP cover, relative shelter volume, juvenile indicator, and
rubble. This is implementation is limited to using outputs that can be provided by
ecological models that only model corals.

# Arguments
- `ltmp_cover` : LTMP Coral Cover with dimensions [timesteps ⋅ locations].
- `relative_shelter_volume` : Relative shelter volume with dimensions [timesteps ⋅ locations].
- `juvenile_indicator`: Juvenile Indicator with dimensions [timesteps ⋅ locations].
- `rubble` : Relative rubble cover with dimensions [timesteps ⋅ locations].
- `out_rci` : Output RCI buffer with dimensions [timesteps ⋅ locations].

# References
1. Ryan F. Heneghan, Gabriela Scheufele, Yves-Marie Bozec et al. A framework to inform
    economic valuation of non-use benefits from coral-reef intervention efforts, 02 October
    2025, PREPRINT (Version 1) available at Research Square
    [https://doi.org/10.21203/rs.3.rs-7644150/v1]
"""
function reef_condition_index!(
    ltmp_cover::AbstractArray{<:AbstractFloat,2},
    relative_shelter_volume::AbstractArray{<:AbstractFloat,2},
    juvenile_indicator::AbstractArray{<:AbstractFloat,2},
    rubble::AbstractArray{<:AbstractFloat,2},
    out_rci::AbstractArray{<:AbstractFloat,2};
)::Nothing
    RCI_CRIT = eltype(ltmp_cover)[
        0.05 0.15 0.25 0.35 0.45;  # relative cover thresholds
        0.15 0.25 0.30 0.35 0.45;  # shelter volume thresholds
        0.15 0.20 0.25 0.30 0.35;  # coral juveniles thresholds
        0.70 0.75 0.80 0.85 0.90  # rubble thresholds
    ]

    CRIT_VAL = eltype(ltmp_cover)[0.1, 0.3, 0.5, 0.7, 0.9]

    CRITERIA_THRESHOLD = 2

    out_rci .= 0.1
    metrics_met = zeros(eltype(ltmp_cover), size(ltmp_cover)...)
    pass_mask = falses(size(ltmp_cover)...)
    for i in 1:size(CRIT_VAL, 1)
        @. metrics_met = (ltmp_cover >= RCI_CRIT[1, i]) +
                         (relative_shelter_volume >= RCI_CRIT[2, i]) +
                         (juvenile_indicator >= RCI_CRIT[3, i]) +
                         ((1.0 .- rubble) >= RCI_CRIT[4, i])
        @. pass_mask = metrics_met >= CRITERIA_THRESHOLD

        out_rci[pass_mask] .= CRIT_VAL[i]
    end


    return nothing
end

"""
    reef_condition_index(relative_cover::AbstractArray{T,2}, shelter_volume::AbstractArray{T,2}, relative_juveniles::AbstractArray{T,2}, rubble::AbstractArray{T,2})::AbstractArray{T,2} where {T<:AbstractFloat}

Calculate the Reef Condition Index (RCI).

This method uses four inputs: LTMP cover, relative shelter volume, juvenile indicator, and rubble.

The RCI is a categorical index assessing the overall health and condition of a reef location
based on five key ecological metrics. The index assigns a discrete score (0.1, 0.3, 0.5, 0.7, or 0.9)
representing categories from "Very Poor" to "Very Good".

For each input there are five levels of condition ranging from very poor to very good. COTS
and Rubble Cover is inverted, where high values indicate worse condition. Then the location
is assigned a condition of very poor to very good if that location meets at least 2 of the
metrics condition criteria. The condition level is then assigned a numerical value based on
its categorisation.

# Arguments
- `ltmp_cover` : Relative Coral Cover with dimensions [timesteps ⋅ locations].
- `relative_shelter_volume` : Relative shelter volume with dimensions [timesteps ⋅ locations].
- `juvenile_indicator`: Relative juvenile cover with dimensions [timesteps ⋅ locations].
- `rubble` : Relative rubble cover with dimensions [timesteps ⋅ locations].

# Returns
Output RCI buffer with dimensions [timesteps ⋅ locations].

# References
1. Ryan F. Heneghan, Gabriela Scheufele, Yves-Marie Bozec et al. A framework to inform
    economic valuation of non-use benefits from coral-reef intervention efforts, 02 October
    2025, PREPRINT (Version 1) available at Research Square
    [https://doi.org/10.21203/rs.3.rs-7644150/v1]
"""
function reef_condition_index(
    relative_cover::AbstractArray{<:AbstractFloat,2},
    shelter_volume::AbstractArray{<:AbstractFloat,2},
    relative_juveniles::AbstractArray{<:AbstractFloat,2},
    rubble::AbstractArray{<:AbstractFloat,2}
)::AbstractArray{<:AbstractFloat,2}
    out_rci::Array{eltype(relative_cover),2} = zeros(
        eltype(relative_cover), size(relative_cover)...
    )
    reef_condition_index!(
        relative_cover, shelter_volume, relative_juveniles, rubble, out_rci
    )

    return out_rci
end

"""
    reef_biodiversity_condition_index!(rc::AbstractArray{T,2}, cd::AbstractArray{T,2}, sv::AbstractArray{T,2}, out_rbci::AbstractArray{T,2})::Nothing where {T<:AbstractFloat}

Calculate the Reef Biodiversity Condition Index (RBCI). This is implementation is limited to
using outputs that can be provided by coral ecology models.

# Arguments
- `rc` : Relative coral cover.
- `cd` : Coral diversity.
- `sv` : Relative shelter volume.
- `out_rbci` : Output array buffer for the RCI.
"""
function reef_biodiversity_condition_index!(
    rc::AbstractArray{<:AbstractFloat,2},
    cd::AbstractArray{<:AbstractFloat,2},
    sv::AbstractArray{<:AbstractFloat,2},
    out_rbci::AbstractArray{<:AbstractFloat,2}
)::Nothing
    out_rbci .= clamp.((rc .+ cd .+ sv) ./ 3.0, 0.0, 1.0)

    return nothing
end

"""
    reef_biodiversity_condition_index(relative_cover::AbstractArray{T,2}, coral_diversity::AbstractArray{T,2}, shelter_volume::AbstractArray{T,2})::AbstractArray{T,2} where {T<:AbstractFloat}

Calculate the Reef Biodiversity Condition Index (RBCI). The RBCI is simply the average of
relative cover (RC), coral diversity (CD), and shelter volume (SV). Given as

```math
\\begin{align*}
\\text{RBCI} = \\frac{\\text{RC} + \\text{CD} + \\text{SV}}{3}.
\\end{align*}
```
This is implementation is limited to using outputs that can be provided by coral ecology
models.

# Arguments
- `relative_cover` : Relative coral cover.
- `coral_diversity` : Coral diversity.
- `shelter_volume` : Relative shelter volume.

# Returns
A 2D array of the Reef Condition Index.
"""
function reef_biodiversity_condition_index(
    relative_cover::AbstractArray{<:AbstractFloat,2},
    coral_diversity::AbstractArray{<:AbstractFloat,2},
    shelter_volume::AbstractArray{<:AbstractFloat,2}
)::AbstractArray{<:AbstractFloat,2}
    rc_size = size(relative_cover)
    if (rc_size != size(coral_diversity)) || (rc_size != size(shelter_volume))
        throw(DimensionMismatch("All input metric arrays must have the same dimensions."))
    end

    out_rbci = zeros(Float64, size(relative_cover))
    reef_biodiversity_condition_index!(
        relative_cover, coral_diversity, shelter_volume, out_rbci
    )

    return out_rbci
end

"""
    reef_tourism_index_no_rubble!(relative_cover::AbstractArray{T,2}, coral_evenness::AbstractArray{T,2}, relative_shelter_volume::AbstractArray{T,2}, relative_juveniles::AbstractArray{T,2}, out_rti::AbstractArray{T,2})::Nothing where {T<:AbstractFloat}

Calculate the Reef Tourism Index (RTI) for a single scenario without rubble as input. This
version of the RTI is for ecological models that do not model rubble.

This method uses four inputs: relative cover, coral evenness, relative shelter volume, and relative juveniles.

# Arguments
- `relative_cover` : Relative coral cover with dimensions [timesteps ⋅ locations].
- `coral_evenness` : Coral evenness with dimensions [timesteps ⋅ locations].
- `relative_shelter_volume` : Relative shelter volume with dimensions [timesteps ⋅ locations].
- `relative_juveniles` : Relative juvenile cover with dimensions [timesteps ⋅ locations].
- `out_rti` : Output array buffer for the RTI with dimensions [timesteps ⋅ locations].
"""
function reef_tourism_index_no_rubble!(
    relative_cover::AbstractArray{<:AbstractFloat,2},
    coral_evenness::AbstractArray{<:AbstractFloat,2},
    relative_shelter_volume::AbstractArray{<:AbstractFloat,2},
    relative_juveniles::AbstractArray{<:AbstractFloat,2},
    out_rti::AbstractArray{<:AbstractFloat,2}
)::Nothing
    # Intercept and coefficient resulting from regressions.
    intcp = 0.47947
    rc_coef = 0.12764
    evenness_coef = 0.31946
    sv_coef = 0.11676
    jv_coef = -0.0036065

    out_rti .= intcp .+
        (rc_coef .* relative_cover) .+
        (evenness_coef .* coral_evenness) .+
        (sv_coef .* relative_shelter_volume) .+
        (jv_coef .* relative_juveniles)

    out_rti .= round.(clamp.(out_rti, 0.1, 0.9); digits=2)

    return nothing
end

"""
    reef_tourism_index_no_rubble(relative_cover::AbstractArray{T,2}, coral_evenness::AbstractArray{T,2}, relative_shelter_volume::AbstractArray{T,2}, relative_juveniles::AbstractArray{T,2})::AbstractArray{T,2} where {T<:AbstractFloat}

Calculate the Reef Tourism Index (RTI) for a single scenario, without rubble as input. This
version of the RTI is for ecological models that do not model rubble.

This method uses four inputs: relative cover, coral evenness, relative shelter volume, and relative juveniles.

The RTI is a continuous version of the Reef Condition Index, fitted with a linear regression model.

# Arguments
- `relative_cover` : Relative coral cover with dimensions [timesteps ⋅ locations].
- `coral_evenness` : Coral evenness with dimensions [timesteps ⋅ locations].
- `relative_shelter_volume` : Relative shelter volume with dimensions [timesteps ⋅ locations].
- `relative_juveniles` : Relative juvenile cover with dimensions [timesteps ⋅ locations].

# Returns
A 2D array of the Reef Tourism Index with dimensions [timesteps ⋅ locations].
"""
function reef_tourism_index_no_rubble(
    relative_cover::AbstractArray{<:AbstractFloat,2},
    coral_evenness::AbstractArray{<:AbstractFloat,2},
    relative_shelter_volume::AbstractArray{<:AbstractFloat,2},
    relative_juveniles::AbstractArray{<:AbstractFloat,2}
)::AbstractArray{<:AbstractFloat,2}
    rc_size = size(relative_cover)
    if (rc_size != size(coral_evenness)) || (rc_size != size(relative_shelter_volume)) || (rc_size != size(relative_juveniles))
        throw(DimensionMismatch("All input metric arrays must have the same dimensions."))
    end

    out_rti = zeros(Float64, rc_size)
    reef_tourism_index_no_rubble!(
        relative_cover,
        coral_evenness,
        relative_shelter_volume,
        relative_juveniles,
        out_rti
    )

    return out_rti
end

"""
    reef_tourism_index!(relative_cover::AbstractArray{T,2}, shelter_volume::AbstractArray{T,2}, relative_juveniles::AbstractArray{T,2}, rubble::AbstractArray{T,2}, out_rti::AbstractArray{T,2})::Nothing where {T<:AbstractFloat}

Calculate the Reef Tourism Index (RTI) for a single scenario.

This method uses four inputs: relative cover, shelter volume, relative juveniles, and rubble.

# Arguments
- `relative_cover` : Relative coral cover with dimensions [timesteps ⋅ locations].
- `shelter_volume` : Relative shelter volume with dimensions [timesteps ⋅ locations].
- `relative_juveniles` : Relative juvenile cover with dimensions [timesteps ⋅ locations].
- `rubble` : Rubble as a proportion of location area with dimensions [timesteps ⋅ locations].
- `out_rti` : Output array buffer for the RTI with dimensions [timesteps ⋅ locations].
"""
function reef_tourism_index!(
    relative_cover::AbstractArray{<:AbstractFloat,2},
    shelter_volume::AbstractArray{<:AbstractFloat,2},
    relative_juveniles::AbstractArray{<:AbstractFloat,2},
    rubble::AbstractArray{<:AbstractFloat,2},
    out_rti::AbstractArray{<:AbstractFloat,2}
)::Nothing
    # Intercept and coefficient resulting from regressions.
    intcp = -0.871
    rc_coef = 0.7678
    sv_coef = 0.2945
    jv_coef = 0.8371
    rubble_coef = 0.7764

    out_rti .= intcp .+
        (rc_coef .* relative_cover) .+
        (sv_coef .* shelter_volume) .+
        (jv_coef .* relative_juveniles) .+
        (rubble_coef .* ( 1 .- rubble))

    out_rti .= round.(clamp.(out_rti, 0.1, 0.9); digits=2)

    return nothing
end

"""
    reef_tourism_index(relative_cover::AbstractArray{T,2}, shelter_volume::AbstractArray{T,2}, relative_juveniles::AbstractArray{T,2}, rubble::AbstractArray{T,2})::AbstractArray{T,2} where {T<:AbstractFloat}

Calculate the Reef Tourism Index (RTI) for a single scenario. The RTI is the Reef
Condition Index made continuous by fitting a linear regression model using relative cover,
shelter volume, relative juveniles, and rubble to underpin it.

This method uses four inputs: relative cover, shelter volume, relative juveniles, and rubble.

# Arguments
- `relative_cover` : Relative Coral Cover with dimensions [timesteps ⋅ locations].
- `shelter_volume` : Relative shelter volume with dimensions [timesteps ⋅ locations].
- `relative_juveniles` : Relative juvenile cover with dimensions [timesteps ⋅ locations].
- `rubble` : Rubble as proportion of location area with dimensions [timesteps ⋅ locations].

# Returns
A 2D array of the Reef Tourism Index with dimensions [timesteps ⋅ locations].
"""
function reef_tourism_index(
    relative_cover::AbstractArray{<:AbstractFloat,2},
    shelter_volume::AbstractArray{<:AbstractFloat,2},
    relative_juveniles::AbstractArray{<:AbstractFloat,2},
    rubble::AbstractArray{<:AbstractFloat,2}
)::AbstractArray{<:AbstractFloat,2}
    rc_size = size(relative_cover)
    if (rc_size != size(shelter_volume)) || (rc_size != size(relative_juveniles)) || (rc_size != size(rubble))
        throw(DimensionMismatch("All input metric arrays must have the same dimensions."))
    end

    out_rti = zeros(eltype(relative_cover), rc_size)
    reef_tourism_index!(
        relative_cover, shelter_volume, relative_juveniles, rubble, out_rti
    )

    return out_rti
end

"""
    reef_fish_index!(rc::AbstractArray, out_rfi::AbstractArray)::Nothing

Calculate the Reef Fish Index (RFI) for a single scenario.

# Arguments
- `rc` : Relative coral cover with dimensions [timesteps ⋅ locations].
- `out_rfi` : Output array buffer for the RFI with dimensions [timesteps ⋅ locations].
"""
function reef_fish_index!(
    rc::AbstractArray{<:AbstractFloat,2},
    out_rfi::AbstractArray{<:AbstractFloat,2}
)::Nothing
    out_rfi .= 0.01 .* (-1623.6 .+ 1883.3 .* (1.232 .+ 0.007476 .* (rc .* 100.0)))
    out_rfi .= round.(out_rfi; digits=2)

    return nothing
end

"""
    reef_fish_index(relative_cover::AbstractArray{T,2},)::AbstractArray{T,2} where {T<:Real}

Calculate the Reef Fish Index (RFI) for a single scenario. The RFI is composed
of two linear regressions mapping relative coral cover to structural complexity and finally,
structural complexity to fish biomass. The index is based off figure 4a and 6b in
Graham et al., 2013 [1]. RFI (kg/km²) as a functional of relative cover (``x``) is given as

```math
\\begin{align*}
\\text{SC}(x) &= 1.232 .+ 0.007476 ⋅ x ⋅ 100\\\\
\\text{RFI} &= 0.01 ⋅ (-1623.6 + 1883.3 ⋅ SC),
\\end{align*}
```
where SC is the structural complexity of the location.

# Arguments
- `relative_cover` : Relative coral cover with dimensions [timesteps ⋅ locations].

# Returns
A 2D array of the Reef Fish Index with dimensions [timesteps ⋅ locations].

# References
1. Graham, N.A.J., Nash, K.L. The importance of structural complexity in coral reef
    ecosystems. Coral Reefs 32, 315–326 (2013). https://doi.org/10.1007/s00338-012-0984-y
"""
function reef_fish_index(
    relative_cover::AbstractArray{<:AbstractFloat,2},
)::AbstractArray{<:AbstractFloat,2}
    out_rfi = zeros(eltype(relative_cover), size(relative_cover))
    reef_fish_index!(relative_cover, out_rfi)

    return out_rfi
end

"    reef_condition_index!(relative_cover::AbstractArray{T,2}, coral_evenness::AbstractArray{T,2}, relative_shelter_volume::AbstractArray{T,2}, relative_juveniles::AbstractArray{T,2}, out_rci::AbstractArray{T,2})::Nothing where {T<:AbstractFloat}

Calculate the Reef Condition Index (RCI).

This method uses four inputs: relative cover, coral evenness, relative shelter volume, and relative juveniles.

# Arguments
- `relative_cover` : Relative Coral Cover with dimensions [timesteps ⋅ locations].
- `relative_shelter_volume` : Relative shelter volume with dimensions [timesteps ⋅ locations].
- `relative_juveniles`: Relative juvenile cover with dimensions [timesteps ⋅ locations].
- `out_rci` : Output RCI buffer with dimensions [timesteps ⋅ locations].
"
function reef_condition_index!(
    relative_cover::AbstractArray{T,2},
    relative_shelter_volume::AbstractArray{T,2},
    relative_juveniles::AbstractArray{T,2},
    out_rci::AbstractArray{T,2};
)::Nothing where {T<:AbstractFloat}
    RCI_CRIT = T[
        0.05 0.15 0.25 0.35 0.45;  # relative cover thresholds
        0.15 0.25 0.30 0.35 0.45;  # shelter volume thresholds
        0.15 0.20 0.25 0.30 0.35   # coral juveniles thresholds
    ]

    CRIT_VAL = T[0.1, 0.3, 0.5, 0.7, 0.9]

    CRITERIA_THRESHOLD = 2

    out_rci .= 0.1
    metrics_met = zeros(T, size(relative_cover)...)
    pass_mask = falses(size(relative_cover)...)
    for i in 1:size(CRIT_VAL, 1)
        @. metrics_met = (relative_cover >= RCI_CRIT[1, i]) +
                         (relative_shelter_volume >= RCI_CRIT[2, i]) +
                         (relative_juveniles >= RCI_CRIT[3, i])
        @. pass_mask = metrics_met >= CRITERIA_THRESHOLD

        out_rci[pass_mask] .= CRIT_VAL[i]
    end


    return nothing
end

"""
    reef_condition_index(relative_cover::AbstractArray{T,2}, coral_evenness::AbstractArray{T,2}, relative_shelter_volume::AbstractArray{T,2}, relative_juveniles::AbstractArray{T,2})::AbstractArray{T,2} where {T<:AbstractFloat}

Calculate the Reef Condition Index (RCI).

This method uses four inputs: relative cover, coral evenness, relative shelter volume, and relative juveniles.

The RCI is a categorical index assessing the overall health and condition of a reef location
based on four key ecological metrics. The index assigns a discrete score (0.1, 0.3, 0.5, 0.7, or 0.9)
representing categories from "Very Poor" to "Very Good".

For each input there are five levels of condition ranging from very poor to very good. Then
the location is assigned a condition of very poor to very good if that location meets 60%
of the metrics condition criteria. The condition level is then assigned a numerical value
based on its categorisation

0.9: Very Good
0.7: Good
0.5: Fair
0.3: Poor
0.1: Very Poor

# Arguments
- `relative_cover` : Relative Coral Cover with dimensions [timesteps ⋅ locations].
- `relative_shelter_volume` : Relative shelter volume with dimensions [timesteps ⋅ locations].
- `relative_juveniles`: Relative juvenile cover with dimensions [timesteps ⋅ locations].

# Returns
Output RCI buffer with dimensions [timesteps ⋅ locations].
"""
function reef_condition_index(
    relative_cover::AbstractArray{T,2},
    relative_shelter_volume::AbstractArray{T,2},
    relative_juveniles::AbstractArray{T,2}
)::AbstractArray{T,2} where {T<:AbstractFloat}
    rc_size = size(relative_cover)
    if (rc_size != size(coral_evenness)) || (rc_size != size(relative_shelter_volume)) || (rc_size != size(relative_juveniles))
        throw(DimensionMismatch("All input metric arrays must have the same dimensions."))
    end

    out_rci = zeros(Float64, size(relative_cover))
    reef_condition_index!(
        relative_cover,
        relative_shelter_volume,
        relative_juveniles,
        out_rci
    )

    return out_rci
end


"""
    reef_condition_index!(relative_cover::AbstractArray{T,2}, relative_shelter_volume::AbstractArray{T,2}, relative_juveniles::AbstractArray{T,2}, cots::AbstractArray{T,2}, rubble::AbstractArray{T,2}, out_rci::AbstractArray{T,2})::Nothing where {T<:AbstractFloat}

Calculate the Reef Condition Index (RCI).

This method uses five inputs: relative cover, relative shelter volume, relative juveniles, COTS abundance, and rubble.

# Arguments
- `relative_cover` : Relative Coral Cover with dimensions [timesteps ⋅ locations].
- `relative_shelter_volume` : Relative shelter volume with dimensions [timesteps ⋅ locations].
- `relative_juveniles`: Relative juvenile cover with dimensions [timesteps ⋅ locations].
- `rubble` : Relative rubble cover with dimensions [timesteps ⋅ locations].
- `out_rci` : Output RCI buffer with dimensions [timesteps ⋅ locations].
"""
function reef_condition_index!(
    relative_cover::AbstractArray{T,2},
    relative_shelter_volume::AbstractArray{T,2},
    relative_juveniles::AbstractArray{T,2},
    rubble::AbstractArray{T,2},
    out_rci::AbstractArray{T,2};
)::Nothing where {T<:AbstractFloat}
    RCI_CRIT = T[
        0.05 0.15 0.25 0.35 0.45;  # relative cover thresholds
        0.15 0.25 0.30 0.35 0.45;  # shelter volume thresholds
        0.15 0.20 0.25 0.30 0.35;  # coral juveniles thresholds
        0.70 0.75 0.80 0.85 0.90  # rubble thresholds
    ]

    CRIT_VAL = T[0.1, 0.3, 0.5, 0.7, 0.9]

    CRITERIA_THRESHOLD = 2

    out_rci .= 0.1
    metrics_met = zeros(T, size(relative_cover)...)
    pass_mask = falses(size(relative_cover)...)
    for i in 1:size(CRIT_VAL, 1)
        @. metrics_met = (relative_cover >= RCI_CRIT[1, i]) +
                         (relative_shelter_volume >= RCI_CRIT[2, i]) +
                         (relative_juveniles >= RCI_CRIT[3, i]) +
                         ((1.0 .- rubble) >= RCI_CRIT[4, i])
        @. pass_mask = metrics_met >= CRITERIA_THRESHOLD

        out_rci[pass_mask] .= CRIT_VAL[i]
    end


    return nothing
end

"""
    reef_condition_index(relative_cover::AbstractArray{T,2}, shelter_volume::AbstractArray{T,2}, relative_juveniles::AbstractArray{T,2}, cots::AbstractArray{T,2}, rubble::AbstractArray{T,2})::AbstractArray{T,2} where {T<:AbstractFloat}

Calculate the Reef Condition Index (RCI).

This method uses five inputs: relative cover, relative shelter volume, relative juveniles, COTS abundance, and rubble.

The RCI is a categorical index assessing the overall health and condition of a reef location
based on five key ecological metrics. The index assigns a discrete score (0.1, 0.3, 0.5, 0.7, or 0.9)
representing categories from "Very Poor" to "Very Good".

For each input there are five levels of condition ranging from very poor to very good. COTS
and Rubble Cover is inverted, where high values indicate worse condition. Then the location
is assigned a condition of very poor to very good if that location meets 60% of the metrics
condition criteria. The condition level is then assigned a numerical value based on its
categorisation

0.9: Very Good
0.7: Good
0.5: Fair
0.3: Poor
0.1: Very Poor

# Arguments
- `relative_cover` : Relative Coral Cover with dimensions [timesteps ⋅ locations].
- `relative_shelter_volume` : Relative shelter volume with dimensions [timesteps ⋅ locations].
- `relative_juveniles`: Relative juvenile cover with dimensions [timesteps ⋅ locations].
- `rubble` : Relative rubble cover with dimensions [timesteps ⋅ locations].

# Returns
Output RCI buffer with dimensions [timesteps ⋅ locations].
"""
function reef_condition_index(
    relative_cover::AbstractArray{T,2},
    shelter_volume::AbstractArray{T,2},
    relative_juveniles::AbstractArray{T,2},
    rubble::AbstractArray{T,2}
)::AbstractArray{T,2} where {T<:AbstractFloat}
    out_rci::Array{T,2} = zeros(T, size(relative_cover)...)
    reef_condition_index!(
        relative_cover, shelter_volume, relative_juveniles, rubble, out_rci
    )

    return out_rci
end

"""
    reef_biodiversity_condition_index!(rc::AbstractArray{T,2}, cd::AbstractArray{T,2}, sv::AbstractArray{T,2}, out_rbci::AbstractArray{T,2})::Nothing where {T<:AbstractFloat}

Calculate the Reef Biodiversity Condition Index (RBCI).

# Arguments
- `rc` : Relative coral cover.
- `cd` : Coral diversity.
- `sv` : Relative shelter volume.
- `out_rbci` : Output array buffer for the RCI.
"""
function reef_biodiversity_condition_index!(
    rc::AbstractArray{T,2},
    cd::AbstractArray{T,2},
    sv::AbstractArray{T,2},
    out_rbci::AbstractArray{T,2}
)::Nothing where {T<:AbstractFloat}
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

# Arguments
- `relative_cover` : Relative coral cover.
- `coral_diversity` : Coral diversity.
- `shelter_volume` : Relative shelter volume.

# Returns
A 2D array of the Reef Condition Index.
"""
function reef_biodiversity_condition_index(
    relative_cover::AbstractArray{T,2},
    coral_diversity::AbstractArray{T,2},
    shelter_volume::AbstractArray{T,2}
)::AbstractArray{T,2} where {T<:AbstractFloat}
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
    reef_tourism_index!(relative_cover::AbstractArray{T,2}, coral_evenness::AbstractArray{T,2}, relative_shelter_volume::AbstractArray{T,2}, relative_juveniles::AbstractArray{T,2}, out_rti::AbstractArray{T,2})::Nothing where {T<:AbstractFloat}

Calculate the Reef Tourism Index (RTI) for a single scenario.

This method uses four inputs: relative cover, coral evenness, relative shelter volume, and relative juveniles.

# Arguments
- `relative_cover` : Relative coral cover with dimensions [timesteps ⋅ locations].
- `coral_evenness` : Coral evenness with dimensions [timesteps ⋅ locations].
- `relative_shelter_volume` : Relative shelter volume with dimensions [timesteps ⋅ locations].
- `relative_juveniles` : Relative juvenile cover with dimensions [timesteps ⋅ locations].
- `out_rti` : Output array buffer for the RTI with dimensions [timesteps ⋅ locations].
"""
function reef_tourism_index!(
    relative_cover::AbstractArray{T,2},
    coral_evenness::AbstractArray{T,2},
    relative_shelter_volume::AbstractArray{T,2},
    relative_juveniles::AbstractArray{T,2},
    out_rti::AbstractArray{T,2}
)::Nothing where {T<:AbstractFloat}
    # Intercept and coefficient resulting from regressions.
    intcp::T = 0.47947
    rc_coef::T = 0.12764
    evenness_coef::T = 0.31946
    sv_coef::T = 0.11676
    jv_coef::T = -0.0036065

    out_rti .= intcp .+
        (rc_coef .* relative_cover) .+
        (evenness_coef .* coral_evenness) .+
        (sv_coef .* relative_shelter_volume) .+
        (jv_coef .* relative_juveniles)

    out_rti .= round.(clamp.(out_rti, 0.1, 0.9); digits=2)

    return nothing
end

"""
    reef_tourism_index(relative_cover::AbstractArray{T,2}, coral_evenness::AbstractArray{T,2}, relative_shelter_volume::AbstractArray{T,2}, relative_juveniles::AbstractArray{T,2})::AbstractArray{T,2} where {T<:AbstractFloat}

Calculate the Reef Tourism Index (RTI) for a single scenario.

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
function reef_tourism_index(
    relative_cover::AbstractArray{T,2},
    coral_evenness::AbstractArray{T,2},
    relative_shelter_volume::AbstractArray{T,2},
    relative_juveniles::AbstractArray{T,2}
)::AbstractArray{T,2} where {T<:AbstractFloat}
    rc_size = size(relative_cover)
    if (rc_size != size(coral_evenness)) || (rc_size != size(relative_shelter_volume)) || (rc_size != size(relative_juveniles))
        throw(DimensionMismatch("All input metric arrays must have the same dimensions."))
    end

    out_rti = zeros(Float64, rc_size)
    reef_tourism_index!(
        relative_cover,
        coral_evenness,
        relative_shelter_volume,
        relative_juveniles,
        out_rti
    )

    return out_rti
end

"""
    reef_tourism_index!(relative_cover::AbstractArray{T,2}, shelter_volume::AbstractArray{T,2}, relative_juveniles::AbstractArray{T,2}, cots::AbstractArray{T,2}, rubble::AbstractArray{T,2}, out_rti::AbstractArray{T,2})::Nothing where {T<:AbstractFloat}

Calculate the Reef Tourism Index (RTI) for a single scenario.

This method uses five inputs: relative cover, shelter volume, relative juveniles, COTS abundance, and rubble.

# Arguments
- `relative_cover` : Relative coral cover with dimensions [timesteps ⋅ locations].
- `shelter_volume` : Relative shelter volume with dimensions [timesteps ⋅ locations].
- `relative_juveniles` : Relative juvenile cover with dimensions [timesteps ⋅ locations].
- `cots` : COTS abundance as a count with dimensions [timesteps ⋅ locations].
- `rubble` : Rubble as a proportion of location area with dimensions [timesteps ⋅ locations].
- `out_rti` : Output array buffer for the RTI with dimensions [timesteps ⋅ locations].
"""
function reef_tourism_index!(
    relative_cover::AbstractArray{T,2},
    shelter_volume::AbstractArray{T,2},
    relative_juveniles::AbstractArray{T,2},
    cots::AbstractArray{T,2},
    rubble::AbstractArray{T,2},
    out_rti::AbstractArray{T,2}
)::Nothing where {T<:AbstractFloat}
    # Intercept and coefficient resulting from regressions.
    intcp::T = 0.47947
    rc_coef::T = 0.7678
    sv_coef::T = 0.2945
    jv_coef::T = 0.8371
    cots_coef::T = 0.2822
    rubble_coef::T = 0.7764

    out_rti .= intcp .+
        (rc_coef .* relative_cover) .+
        (sv_coef .* shelter_volume) .+
        (jv_coef .* relative_juveniles) .+
        (cots_coef .* cots) .+
        (rubble_coef .* rubble)

    out_rti .= round.(clamp.(out_rti, 0.1, 0.9); digits=2)

    return nothing
end

"""
    reef_tourism_index(relative_cover::AbstractArray{T,2}, shelter_volume::AbstractArray{T,2}, relative_juveniles::AbstractArray{T,2}, cots::AbstractArray{T,2}, rubble::AbstractArray{T,2})::AbstractArray{T,2} where {T<:AbstractFloat}

Calculate the Reef Tourism Index (RTI) for a single scenario. The RTI is the Reef
Condition Index made continuous by fitting a linear regression model using relative cover,
shelter volume, relative juveniles, cots abundance and rubble to underpin it.

This method uses five inputs: relative cover, shelter volume, relative juveniles, COTS abundance, and rubble.

# Arguments
- `relative_cover` : Relative Coral Cover with dimensions [timesteps ⋅ locations].
- `shelter_volume` : Relative shelter volume with dimensions [timesteps ⋅ locations].
- `relative_juveniles` : Relative juvenile cover with dimensions [timesteps ⋅ locations].
- `cots` : COTS abundance with dimensions [timesteps ⋅ locations].
- `rubble` : Rubble as proportion of location area with dimensions [timesteps ⋅ locations].

# Returns
A 2D array of the Reef Tourism Index with dimensions [timesteps ⋅ locations].
"""
function reef_tourism_index(
    relative_cover::AbstractArray{T,2},
    shelter_volume::AbstractArray{T,2},
    relative_juveniles::AbstractArray{T,2},
    cots::AbstractArray{T,2},
    rubble::AbstractArray{T,2}
)::AbstractArray{T,2} where {T<:AbstractFloat}
    rc_size = size(relative_cover)
    if (rc_size != size(shelter_volume)) || (rc_size != size(relative_juveniles)) || (rc_size != size(cots)) || (rc_size != size(cots)) || (rc_size != size(rubble))
        throw(DimensionMismatch("All input metric arrays must have the same dimensions."))
    end

    out_rti = zeros(Float64, rc_size)
    reef_tourism_index!(
        relative_cover, shelter_volume, relative_juveniles, cots, rubble, out_rti
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
    rc::AbstractArray{T,2},
    out_rfi::AbstractArray{T,2}
)::Nothing where {T<:Real}
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
    relative_cover::AbstractArray{T,2},
)::AbstractArray{T,2} where {T<:Real}
    out_rfi = zeros(Float64, size(relative_cover))
    reef_fish_index!(relative_cover, out_rfi)

    return out_rfi
end

"""
    _reef_biodiversity_condition_index!(rc::Array{T,2}, cd::AbstractArray{T,2}, sv::AbstractArray{T,2}, out_rci::AbstractArray{T,2})::Nothing where {T<:AbstractFloat}

Calculate the Reef Biodiversity Condition Index (RBCI).

# Arguments
- `rc` : Relative coral cover.
- `cd` : Coral diversity.
- `sv` : Relative shelter volume.
- `out_rci` : Output array buffer for the RCI.
"""
function _reef_biodiversity_condition_index!(
    rc::Array{T,2},
    cd::AbstractArray{T,2},
    sv::AbstractArray{T,2},
    out_rci::AbstractArray{T,2}
)::Nothing where {T<:AbstractFloat}
    out_rci .= clamp.((rc .+ cd .+ sv) ./ 3.0, 0.0, 1.0)

    return nothing
end

"""
    reef_biodiversity_condition_index(relative_cover::Array{T,2}, coral_diversity::Array{T,2}, shelter_volume::Array{T,2})::Array{T,2} where {T<:AbstractFloat}

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
    relative_cover::Array{T,2},
    coral_diversity::Array{T,2},
    shelter_volume::Array{T,2}
)::Array{T,2} where {T<:AbstractFloat}
    rc_size = size(relative_cover)
    if (rc_size != size(coral_diversity)) || (rc_size != size(shelter_volume))
        throw(DimensionMismatch("All input metric arrays must have the same dimensions."))
    end

    out_rci = zeros(Float64, size(relative_cover))
    _reef_biodiversity_condition_index!(
        relative_cover, coral_diversity, shelter_volume, out_rci
    )

    return out_rci
end

"""
    _reef_tourism_index!(relative_cover::Array{T,2}, shelter_volume::Array{T,2}, relative_juveniles::Array{T,2}, cots::Array{T,2}, rubble::Array{T,2}, out_rti::Array{T,2})::Nothing where {T<:AbstractFloat}

Calculate the Reef Tourism Index (RTI) for a single scenario.

# Arguments
- `relative_cover` : Relative coral cover.
- `shelter_volume` : Relative shelter volume.
- `relative_juveniles` : Relative juvenile cover.
- `cots` : COTS abundance.
- `rubble` : Rubble.
- `out_rti` : Output array buffer for the RTI.
"""
function _reef_tourism_index!(
    relative_cover::Array{T,2},
    shelter_volume::Array{T,2},
    relative_juveniles::Array{T,2},
    cots::Array{T,2},
    rubble::Array{T,2},
    out_rti::Array{T,2}
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
    reef_tourism_index(relative_cover::Array{T,2}, shelter_volume::Array{T,2}, relative_juveniles::Array{T,2}, cots::Array{T,2}, rubble::Array{T,2})::Array{T,2} where {T<:AbstractFloat}

Calculate the Reef Tourism Index (RTI) for a single scenario. The RTI is the Reef
Condition Index made continuous by fitting a linear regression model using relative cover,
shelter volume, relative juveniles, cots abundance and rubble to underpin it.

# Arguments
- `relative_cover` : Relative coral cover.
- `shelter_volume` : Relative shelter volume.
- `relative_juveniles` : Relative juvenile cover.
- `cots` : COTS abundance.
- `rubble` : Rubble.

# Returns
A 2D array of the Reef Tourism Index.
"""
function reef_tourism_index(
    relative_cover::Array{T,2},
    shelter_volume::Array{T,2},
    relative_juveniles::Array{T,2},
    cots::Array{T,2},
    rubble::Array{T,2}
)::Array{T,2} where {T<:AbstractFloat}
    rc_size = size(relative_cover)
    if (rc_size != size(shelter_volume)) || (rc_size != size(relative_juveniles)) || (rc_size != size(cots)) || (rc_size != size(cots)) || (rc_size != size(rubble))
        throw(DimensionMismatch("All input metric arrays must have the same dimensions."))
    end

    out_rti = zeros(Float64, rc_size)
    _reef_tourism_index!(
        relative_cover, shelter_volume, relative_juveniles, cots, rubble, out_rti
    )

    return out_rti
end

"""
    _reef_fish_index!(rc::AbstractArray, out_rfi::AbstractArray)::Nothing

Calculate the Reef Fish Index (RFI) for a single scenario.

# Arguments
- `rc` : Relative coral cover.
- `out_rfi` : Output array buffer for the RFI.
"""
function _reef_fish_index!(
    rc::AbstractArray,
    out_rfi::AbstractArray
)::Nothing
    out_rfi .= 0.01 .* (-1623.6 .+ 1883.3 .* (1.232 .+ 0.007476 .* (rc .* 100.0)))
    out_rfi .= round.(out_rfi; digits=2)

    return nothing
end

"""
    reef_fish_index(relative_cover::Array{T,2},)::Array{T,2} where {T<:Real}

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
A 2D array of the Reef Fish Index.

# References
1. Graham, N.A.J., Nash, K.L. The importance of structural complexity in coral reef
    ecosystems. Coral Reefs 32, 315–326 (2013). https://doi.org/10.1007/s00338-012-0984-y
"""
function reef_fish_index(
    relative_cover::Array{T,2},
)::Array{T,2} where {T<:Real}
    out_rfi = zeros(Float64, size(relative_cover))
    _reef_fish_index!(relative_cover, out_rfi)

    return out_rfi
end

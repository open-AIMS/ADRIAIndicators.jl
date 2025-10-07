"""
    relative_cover!(relative_cover::AbstractArray{T,4}, out_relative_cover::AbstractArray{T,2})::Nothing where {T<:Real}

Calculate the relative cover per location by summing over groups and size classes.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `out_relative_cover` : Output array buffer for relative cover with dimensions [timesteps ⋅ locations].
"""
function relative_cover!(
    relative_cover::AbstractArray{T,4},
    out_relative_cover::AbstractArray{T,2}
)::Nothing where {T<:Real}
    # Sum over groups and sizes
    out_relative_cover .= dropdims(sum(relative_cover; dims=(2, 3)); dims=(2, 3))

    return nothing
end

"""
    relative_cover(relative_cover::AbstractArray{T,4})::AbstractArray{T,2} where {T<:Real}

Calculate the relative cover per location by summing over groups and size classes.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].

# Returns
A 2D array of relative cover with dimensions [timesteps ⋅ locations].
"""
function relative_cover(relative_cover::AbstractArray{T,4})::AbstractArray{T,2} where {T<:Real}
    n_timesteps, _, _, n_locations = size(relative_cover)
    out_relative_cover = zeros(T, n_timesteps, n_locations)
    relative_cover!(relative_cover, out_relative_cover)

    return out_relative_cover
end

"""
    relative_taxa_cover!(relative_cover::AbstractArray{T,4}, location_area::AbstractVector{T}, out_relative_taxa_cover::AbstractArray{T,2})::Nothing where {T<:Real}

Calculate the relative taxa cover, summed up across all locations.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `location_area` : The coral habitable area for each location.
- `out_relative_taxa_cover` : Output array buffer with dimensions [timesteps ⋅ groups].
"""
function relative_taxa_cover!(
    relative_cover::AbstractArray{T,4},
    location_area::AbstractVector{T},
    out_relative_taxa_cover::AbstractArray{T,2}
)::Nothing where {T<:Real}
    # Sum over sizes
    group_cover = dropdims(sum(relative_cover; dims=3); dims=3)  # [timesteps, groups, locations]
    absolute_group_cover = group_cover .* reshape(location_area, (1, 1, :))
    total_location_area = sum(location_area)
    out_relative_taxa_cover .= 
        dropdims(sum(absolute_group_cover; dims=3); dims=3) ./ total_location_area

    return nothing
end
"""
    relative_taxa_cover(relative_cover::AbstractArray{T,4}, location_area::AbstractVector{T})::AbstractArray{T,2} where {T<:Real}

Calculate the relative taxa cover, summed up across all locations.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `location_area` : The coral habitable area for each location.

# Returns
A 2D array of relative taxa cover with dimensions [timesteps ⋅ groups].
"""
function relative_taxa_cover(
    relative_cover::AbstractArray{T,4},
    location_area::AbstractVector{T}
)::AbstractArray{T,2} where {T<:Real}
    n_timesteps, n_groups, _, n_locations = size(relative_cover)
    if length(location_area) != n_locations
        throw(
            DimensionMismatch(
                "The number of locations in relative_cover and location_area must match."
            )
        )
    end
    out_relative_taxa_cover = zeros(T, n_timesteps, n_groups)
    relative_taxa_cover!(relative_cover, location_area, out_relative_taxa_cover)

    return out_relative_taxa_cover
end

"""
    relative_loc_taxa_cover!(relative_cover::AbstractArray{T,4}, out_relative_loc_taxa_cover::AbstractArray{T,3})::Nothing where {T<:Real}

Calculate the relative taxa cover for each location.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `out_relative_loc_taxa_cover` : Output array buffer with dimensions [timesteps ⋅ groups ⋅ locations].
"""
function relative_loc_taxa_cover!(
    relative_cover::AbstractArray{T,4},
    out_relative_loc_taxa_cover::AbstractArray{T,3}
)::Nothing where {T<:Real}
    out_relative_loc_taxa_cover .= dropdims(sum(relative_cover; dims=3); dims=3)

    return nothing
end

"""
    relative_loc_taxa_cover(relative_cover::AbstractArray{T,4})::AbstractArray{T,3} where {T<:Real}

Calculate the relative taxa cover for each location.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].

# Returns
A 3D array of relative taxa cover with dimensions [timesteps ⋅ groups ⋅ locations].
"""
function relative_loc_taxa_cover(
    relative_cover::AbstractArray{T,4}
)::AbstractArray{T,3} where {T<:Real}
    n_timesteps, n_groups, _, n_locations = size(relative_cover)
    out_relative_loc_taxa_cover = zeros(T, n_timesteps, n_groups, n_locations)
    relative_loc_taxa_cover!(relative_cover, out_relative_loc_taxa_cover)

    return out_relative_loc_taxa_cover
end

"""
    ltmp_cover!(relative_cover_input::AbstractArray{T,4}, habitable_area::AbstractVector{T}, reef_area::AbstractVector{T}, out_ltmp_cover::AbstractArray{T,2})::Nothing where {T<:Real}

Calculate the LTMP cover for each location. LTMP cover is the proportion of the reef area
occupied by coral relative to the area of the reef, which may include non-habitable area.
The `relative_cover_input` is relative to habitable area.

# Arguments
- `relative_cover_input` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `habitable_area` : The habitable area for each location in m².
- `reef_area` : The total area for each location in m².
- `out_ltmp_cover` : Output array buffer for LTMP cover with dimensions [timesteps ⋅ locations].
"""
function ltmp_cover!(
    relative_cover_input::AbstractArray{T,4},
    habitable_area::AbstractVector{T},
    reef_area::AbstractVector{T},
    out_ltmp_cover::AbstractArray{T,2}
)::Nothing where {T<:Real}
    # First, calculate cover relative to habitable area by summing over groups and sizes
    rel_cover = zeros(T, size(out_ltmp_cover)...)
    relative_cover!(relative_cover_input, rel_cover)

    # Convert relative cover to LTMP cover
    # Note: The location dimension for the 2D rel_cover array is 2
    relative_cover_to_ltmp_cover!(rel_cover, habitable_area, reef_area, 2, out_ltmp_cover)

    return nothing
end

"""
    ltmp_cover(relative_cover_input::AbstractArray{T,4}, habitable_area::AbstractVector{T}, reef_area::AbstractVector{T})::AbstractArray{T,2} where {T<:Real}

Calculate the LTMP cover for each location. LTMP cover is the proportion of the reef area
occupied by coral relative to the area of the reef, which may include non-habitable area.
The `relative_cover_input` is relative to habitable area.

# Arguments
- `relative_cover_input` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `habitable_area` : The habitable area for each location in m².
- `reef_area` : The total area for each location in m².

# Returns
A 2D array of LTMP cover with dimensions [timesteps ⋅ locations].
"""
function ltmp_cover(
    relative_cover_input::AbstractArray{T,4},
    habitable_area::AbstractVector{T},
    reef_area::AbstractVector{T}
)::AbstractArray{T,2} where {T<:Real}
    n_timesteps, _, _, n_locations = size(relative_cover_input)
    if length(habitable_area) != n_locations || length(reef_area) != n_locations
        throw(
            DimensionMismatch(
                "The number of locations in relative_cover_input, habitable_area, and reef_area must match."
            )
        )
    end

    out_ltmp_cover = zeros(T, n_timesteps, n_locations)
    ltmp_cover!(relative_cover_input, habitable_area, reef_area, out_ltmp_cover)

    return out_ltmp_cover
end

"""
    ltmp_taxa_cover!(relative_cover::AbstractArray{T,4}, habitable_area::AbstractVector{T}, reef_area::AbstractVector{T}, out_ltmp_taxa_cover::AbstractArray{T,2})::Nothing where {T<:Real}

Calculate the LTMP taxa cover, summed up across all locations.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `habitable_area` : The habitable area for each location in m².
- `reef_area` : The total area for each location in m².
- `out_ltmp_taxa_cover` : Output array buffer with dimensions [timesteps ⋅ groups].
"""
function ltmp_taxa_cover!(
    relative_cover::AbstractArray{T,4},
    habitable_area::AbstractVector{T},
    reef_area::AbstractVector{T},
    out_ltmp_taxa_cover::AbstractArray{T,2}
)::Nothing where {T<:Real}
    # Get group cover relative to habitable area for each location
    group_cover = dropdims(sum(relative_cover; dims=3); dims=3)  # [timesteps, groups, locations]

    # Convert to absolute cover in m²
    absolute_group_cover = group_cover .* reshape(habitable_area, (1, 1, :))

    # Sum absolute cover over all locations for each group
    total_absolute_group_cover = dropdims(sum(absolute_group_cover; dims=3); dims=3) # [timesteps, groups]

    # Get total reef area
    total_reef_area = sum(reef_area)

    # Divide by total reef area to get LTMP taxa cover
    out_ltmp_taxa_cover .= total_absolute_group_cover ./ total_reef_area

    return nothing
end

"""
    ltmp_taxa_cover(relative_cover::AbstractArray{T,4}, habitable_area::AbstractVector{T}, reef_area::AbstractVector{T})::AbstractArray{T,2} where {T<:Real}

Calculate the LTMP taxa cover, summed up across all locations.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `habitable_area` : The habitable area for each location in m².
- `reef_area` : The total area for each location in m².

# Returns
A 2D array of LTMP taxa cover with dimensions [timesteps ⋅ groups].
"""
function ltmp_taxa_cover(
    relative_cover::AbstractArray{T,4},
    habitable_area::AbstractVector{T},
    reef_area::AbstractVector{T}
)::AbstractArray{T,2} where {T<:Real}
    n_timesteps, n_groups, _, n_locations = size(relative_cover)
    if length(habitable_area) != n_locations || length(reef_area) != n_locations
        throw(
            DimensionMismatch(
                "The number of locations in relative_cover, habitable_area, and reef_area must match."
            )
        )
    end
    out_ltmp_taxa_cover = zeros(T, n_timesteps, n_groups)
    ltmp_taxa_cover!(relative_cover, habitable_area, reef_area, out_ltmp_taxa_cover)

    return out_ltmp_taxa_cover
end

"""
    ltmp_loc_taxa_cover!(relative_cover::AbstractArray{T,4}, habitable_area::AbstractVector{T}, reef_area::AbstractVector{T}, out_ltmp_loc_taxa_cover::AbstractArray{T,3})::Nothing where {T<:Real}

Calculate the LTMP taxa cover for each location.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `habitable_area` : The habitable area for each location in m².
- `reef_area` : The total area for each location in m².
- `out_ltmp_loc_taxa_cover` : Output array buffer with dimensions [timesteps ⋅ groups ⋅ locations].
"""
function ltmp_loc_taxa_cover!(
    relative_cover::AbstractArray{T,4},
    habitable_area::AbstractVector{T},
    reef_area::AbstractVector{T},
    out_ltmp_loc_taxa_cover::AbstractArray{T,3}
)::Nothing where {T<:Real}
    # Get taxa cover relative to habitable area
    relative_loc_taxa_cover!(relative_cover, out_ltmp_loc_taxa_cover)

    # Convert to be relative to reef area for each location
    area_coefficient = reshape(habitable_area ./ reef_area, (1, 1, :))
    out_ltmp_loc_taxa_cover .*= area_coefficient

    return nothing
end

"""
    ltmp_loc_taxa_cover(relative_cover::AbstractArray{T,4}, habitable_area::AbstractVector{T}, reef_area::AbstractVector{T})::AbstractArray{T,3} where {T<:Real}

Calculate the LTMP taxa cover for each location.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `habitable_area` : The habitable area for each location in m².
- `reef_area` : The total area for each location in m².

# Returns
A 3D array of LTMP taxa cover with dimensions [timesteps ⋅ groups ⋅ locations].
"""
function ltmp_loc_taxa_cover(
    relative_cover::AbstractArray{T,4},
    habitable_area::AbstractVector{T},
    reef_area::AbstractVector{T}
)::AbstractArray{T,3} where {T<:Real}
    n_timesteps, n_groups, _, n_locations = size(relative_cover)
    if length(habitable_area) != n_locations || length(reef_area) != n_locations
        throw(
            DimensionMismatch(
                "The number of locations in relative_cover, habitable_area, and reef_area must match."
            )
        )
    end
    out_ltmp_loc_taxa_cover = zeros(T, n_timesteps, n_groups, n_locations)
    ltmp_loc_taxa_cover!(relative_cover, habitable_area, reef_area, out_ltmp_loc_taxa_cover)

    return out_ltmp_loc_taxa_cover
end
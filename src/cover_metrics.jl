"""
    relative_cover!(relative_cover::AbstractArray{T,4}, out_relative_cover::AbstractArray{T,2})::Nothing where {T<:Real}

Calculate the relative cover per location by summing over groups and size classes.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations], relative to habitable area.
- `out_relative_cover` : Output array buffer for relative cover with dimensions [timesteps ⋅ locations].
"""
function relative_cover!(
    relative_cover::AbstractArray{T,4},
    out_relative_cover::AbstractArray{T,2}
)::Nothing where {T<:Real}
    # Sum over groups and sizes (dimensions 2 and 3)
    n_timesteps, _, _, n_locations = size(relative_cover)
    sum!(reshape(out_relative_cover, (n_timesteps, 1, 1, n_locations)), relative_cover)

    return nothing
end

"""
    relative_cover(relative_cover::AbstractArray{T,4})::AbstractArray{T,2} where {T<:Real}

Calculate the relative cover per location by summing over groups and size classes.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations], relative to habitable area.

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
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations], relative to habitable area.
- `location_area` : The coral habitable area for each location.
- `out_relative_taxa_cover` : Output array buffer with dimensions [timesteps ⋅ groups].
"""
function relative_taxa_cover!(
    relative_cover::AbstractArray{T,4},
    location_area::AbstractVector{T},
    out_relative_taxa_cover::AbstractArray{T,2}
)::Nothing where {T<:Real}
    n_timesteps, n_groups, n_sizes, n_locations = size(relative_cover)
    total_area = sum(location_area)

    fill!(out_relative_taxa_cover, zero(T))

    for l in 1:n_locations
        area = location_area[l]
        for s in 1:n_sizes
            for g in 1:n_groups
                for t in 1:n_timesteps
                    out_relative_taxa_cover[t, g] += relative_cover[t, g, s, l] * area
                end
            end
        end
    end

    out_relative_taxa_cover ./= total_area

    return nothing
end
"""
    relative_taxa_cover(relative_cover::AbstractArray{T,4}, location_area::AbstractVector{T})::AbstractArray{T,2} where {T<:Real}

Calculate the relative taxa cover, summed up across all locations.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations], relative to habitable area.
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
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations], relative to habitable area.
- `out_relative_loc_taxa_cover` : Output array buffer with dimensions [timesteps ⋅ groups ⋅ locations].
"""
function relative_loc_taxa_cover!(
    relative_cover::AbstractArray{T,4},
    out_relative_loc_taxa_cover::AbstractArray{T,3}
)::Nothing where {T<:Real}
    n_timesteps, n_groups, _, n_locations = size(relative_cover)
    sum!(reshape(out_relative_loc_taxa_cover, (n_timesteps, n_groups, 1, n_locations)), relative_cover)

    return nothing
end

"""
    relative_loc_taxa_cover(relative_cover::AbstractArray{T,4})::AbstractArray{T,3} where {T<:Real}

Calculate the relative taxa cover for each location.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations], relative to habitable area.

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
    ltmp_cover!(relative_cover::AbstractArray{T,4}, habitable_area::AbstractVector{T}, reef_area::AbstractVector{T}, out_ltmp_cover::AbstractArray{T,2})::Nothing where {T<:Real}

Calculate the LTMP cover for each location. LTMP cover is the proportion of the reef area
occupied by coral relative to the area of the reef, which may include non-habitable area.
The `relative_cover` input is relative to habitable area.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `habitable_area` : The habitable area for each location in m².
- `reef_area` : The total area for each location in m².
- `out_ltmp_cover` : Output array buffer for LTMP cover with dimensions [timesteps ⋅ locations].
"""
function ltmp_cover!(
    relative_cover::AbstractArray{T,4},
    habitable_area::AbstractVector{T},
    reef_area::AbstractVector{T},
    out_ltmp_cover::AbstractArray{T,2}
)::Nothing where {T<:Real}
    # Calculate cover relative to habitable area by summing over groups and sizes
    # directly into out_ltmp_cover
    n_timesteps, _, _, n_locations = size(relative_cover)
    sum!(reshape(out_ltmp_cover, (n_timesteps, 1, 1, n_locations)), relative_cover)

    # Convert relative cover to LTMP cover in-place
    # LTMP = RC * (Habitable Area / Reef Area)
    # The location dimension is the 2nd dimension of out_ltmp_cover
    area_coeff = reshape(habitable_area ./ reef_area, (1, :))
    out_ltmp_cover .*= area_coeff

    return nothing
end

"""
    ltmp_cover(relative_cover::AbstractArray{T,4}, habitable_area::AbstractVector{T}, reef_area::AbstractVector{T})::AbstractArray{T,2} where {T<:Real}

Calculate the LTMP cover for each location. LTMP cover is the proportion of the reef area
occupied by coral relative to the area of the reef, which may include non-habitable area.
The `relative_cover` input is relative to habitable area.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `habitable_area` : The habitable area for each location in m².
- `reef_area` : The total area for each location in m².

# Returns
A 2D array of LTMP cover with dimensions [timesteps ⋅ locations].
"""
function ltmp_cover(
    relative_cover::AbstractArray{T,4},
    habitable_area::AbstractVector{T},
    reef_area::AbstractVector{T}
)::AbstractArray{T,2} where {T<:Real}
    n_timesteps, _, _, n_locations = size(relative_cover)
    if length(habitable_area) != n_locations || length(reef_area) != n_locations
        throw(
            DimensionMismatch(
                "The number of locations in relative_cover, habitable_area, and reef_area must match."
            )
        )
    end

    out_ltmp_cover = zeros(T, n_timesteps, n_locations)
    ltmp_cover!(relative_cover, habitable_area, reef_area, out_ltmp_cover)

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
    n_timesteps, n_groups, n_sizes, n_locations = size(relative_cover)
    total_reef_area = sum(reef_area)

    fill!(out_ltmp_taxa_cover, zero(T))

    for l in 1:n_locations
        h_area = habitable_area[l]
        for s in 1:n_sizes
            for g in 1:n_groups
                for t in 1:n_timesteps
                    out_ltmp_taxa_cover[t, g] += relative_cover[t, g, s, l] * h_area
                end
            end
        end
    end

    out_ltmp_taxa_cover ./= total_reef_area

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
    n_timesteps, n_groups, _, n_locations = size(relative_cover)
    sum!(reshape(out_ltmp_loc_taxa_cover, (n_timesteps, n_groups, 1, n_locations)), relative_cover)

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

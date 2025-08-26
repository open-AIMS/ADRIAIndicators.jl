"""
    relative_loc_cover!(relative_cover::AbstractArray{T,4}, out_relative_cover::AbstractArray{T,2})::Nothing where {T<:Real}

Calculate the relative cover per location by summing over groups and size classes.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `out_relative_cover` : Output array buffer for relative cover with dimensions [timesteps ⋅ locations].
"""
function relative_loc_cover!(
    relative_cover::AbstractArray{T,4},
    out_relative_cover::AbstractArray{T,2}
)::Nothing where {T<:Real}
    # Sum over groups and sizes
    out_relative_cover .= dropdims(sum(relative_cover; dims=(2, 3)); dims=(2, 3))

    return nothing
end

"""
    relative_cover(relative_cover::AbstractArray{T,4})::AbstractArray{T,2} where {T<:Real}

Calculate the relative cover by summing over groups and size classes.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].

# Returns
A 2D array of relative cover with dimensions [timesteps ⋅ locations].
"""
function relative_loc_cover(relative_cover::AbstractArray{T,4})::AbstractArray{T,2} where {T<:Real}
    n_timesteps, _, _, n_locations = size(relative_cover)
    out_relative_cover = zeros(T, n_timesteps, n_locations)
    relative_loc_cover!(relative_cover, out_relative_cover)

    return out_relative_cover
end

"""
    relative_taxa_cover!(relative_cover::AbstractArray{T,4}, k_area::AbstractVector{T}, out_relative_taxa_cover::AbstractArray{T,2})::Nothing where {T<:Real}

Calculate the relative taxa cover, summed up across all locations.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `k_area` : The coral habitable area for each location.
- `out_relative_taxa_cover` : Output array buffer with dimensions [timesteps ⋅ groups].
"""
function relative_taxa_cover!(
    relative_cover::AbstractArray{T,4},
    k_area::AbstractVector{T},
    out_relative_taxa_cover::AbstractArray{T,2}
)::Nothing where {T<:Real}
    # Sum over sizes
    group_cover = dropdims(sum(relative_cover; dims=3); dims=3)  # [timesteps, groups, locations]
    absolute_group_cover = group_cover .* k_area'
    total_k_area = sum(k_area)
    out_relative_taxa_cover .=
        dropdims(sum(absolute_group_cover; dims=3); dims=3) ./ total_k_area

    return nothing
end

"""
    relative_taxa_cover(relative_cover::AbstractArray{T,4}, k_area::AbstractVector{T})::AbstractArray{T,2} where {T<:Real}

Calculate the relative taxa cover, summed up across all locations.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `k_area` : The coral habitable area for each location.

# Returns
A 2D array of relative taxa cover with dimensions [timesteps ⋅ groups].
"""
function relative_taxa_cover(
    relative_cover::AbstractArray{T,4},
    k_area::AbstractVector{T}
)::AbstractArray{T,2} where {T<:Real}
    n_timesteps, n_groups, _, n_locations = size(relative_cover)
    if length(k_area) != n_locations
        throw(
            DimensionMismatch(
                "The number of locations in relative_cover and k_area must match."
            )
        )
    end
    out_relative_taxa_cover = zeros(T, n_timesteps, n_groups)
    relative_taxa_cover!(relative_cover, k_area, out_relative_taxa_cover)

    return out_relative_taxa_cover
end

"""
    relative_cover!(relative_cover::AbstractArray{T,4}, location_area::AbstractVector{T}, out_relative_cover::AbstractVector{T})::Nothing where {T<:AbstractFloat}

For each timestep, calculate the proportion of the entire study area covered by coral.

# Arguments
- `relative_cover` : Raw relative coral cover of dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations]
- `location_areas` : Location area for each location.
- `out_relative_cover` : Vector of dimensions [timesteps].
"""
function relative_cover!(
    relative_cover::AbstractArray{T,4},
    location_area::AbstractVector{T},
    out_relative_cover::AbstractVector{T}
)::Nothing where {T<:AbstractFloat}
    total_area::T = sum(location_area)
    out_relative_cover .=
        dropdims(
            sum(
                relative_cover .* reshape(location_area, (1, 1, 1, -1)); dims=(2, 3, 4)
            ); dims=(2, 3, 4)) ./ total_area

    return nothing
end

"""
    relative_cover(relative_cover::AbstractArray{T,4}, location_area::AbstractVector{T})::AbstractVector{T} where {T<:AbstractFloat}

For each timestep, calculate the proportion of the entire study area covered by coral.

# Arguments
- `relative_cover` : Raw relative coral cover of dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations]
- `location_areas` : Location area for each location.

# Returns
Vector of dimensions [timesteps].
"""
function relative_cover(
    relative_cover::AbstractArray{T,4},
    location_area::AbstractVector{T}
)::AbstractVector{T} where {T<:AbstractFloat}
    n_timesteps = size(relative_cover, 4)
    out_relative_cover = zeros(T, n_timesteps)
    relative_cover!(relative_cover, location_area, out_relative_cover)

    return out_relative_cover
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

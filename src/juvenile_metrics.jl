"""
    relative_juveniles!(relative_cover::AbstractArray{T,4}, is_juvenile::AbstractVector{Bool}, out_relative_juveniles::AbstractArray{T,2})::Nothing where {T<:Real}

Calculate the relative coral cover composed of juveniles.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `is_juvenile` : A boolean vector indicating which size classes are juvenile.
- `out_relative_juveniles` : Output array buffer with dimensions [timesteps ⋅ locations].
"""
function relative_juveniles!(
    relative_cover::AbstractArray{T,4},
    is_juvenile::AbstractVector{Bool},
    out_relative_juveniles::AbstractArray{T,2}
)::Nothing where {T<:Real}
    juvenile_cover = relative_cover[:, :, is_juvenile, :]
    out_relative_juveniles .= dropdims(sum(juvenile_cover; dims=(2, 3)); dims=(2, 3))

    return nothing
end

"""
    relative_juveniles(relative_cover::AbstractArray{T,4}, is_juvenile::AbstractVector{Bool})::AbstractArray{T,2} where {T<:Real}

Calculate the relative coral cover composed of juveniles.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `is_juvenile` : A boolean vector indicating which size classes are juvenile.

# Returns
A 2D array of relative juvenile cover with dimensions [timesteps ⋅ locations].
"""
function relative_juveniles(
    relative_cover::AbstractArray{T,4},
    is_juvenile::AbstractVector{Bool}
)::AbstractArray{T,2} where {T<:Real}
    n_timesteps, _, n_sizes, n_locations = size(relative_cover)
    if length(is_juvenile) != n_sizes
        throw(
            DimensionMismatch(
                "The length of is_juvenile must match the number of size classes in relative_cover."
            )
        )
    end
    out_relative_juveniles = zeros(T, n_timesteps, n_locations)
    relative_juveniles!(relative_cover, is_juvenile, out_relative_juveniles)

    return out_relative_juveniles
end

"""
    relative_taxa_juveniles!(relative_cover::AbstractArray{T,4}, is_juvenile::AbstractArray{Bool}, location_area::AbstractVector{T}, out_relative_taxa_juveniles::AbstractArray{T,2})::Nothing where {T<:AbstractFloat}

Calculate the relative coral cover composed of juveniles over time and functional group.
Returns the output into a preallocated buffer.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `is_juvenile` : A boolean vector indicating which size classes are juvenile.
- `location_area` : Vector containing the area of each location with dimensions [locations].
- `out_relative_taxa_juveniles` : Output array buffer with dimensions [timesteps ⋅ groups].
"""
function relative_taxa_juveniles!(
    relative_cover::AbstractArray{T,4},
    is_juvenile::AbstractArray{Bool},
    location_area::AbstractVector{T},
    out_relative_taxa_juveniles::AbstractArray{T,2}
)::Nothing where {T<:AbstractFloat}
    _is_juveniles = reshape(is_juvenile, (1, 1, :, 1))
    _location_area = reshape(location_area, (1, 1, 1, :))
    out_relative_taxa_juveniles .=
        dropdims(sum(
                relative_cover .* _location_area .* _is_juveniles; dims=(3, 4)
            ); dims=(3, 4)) ./ sum(location_area)

    return nothing
end

"""
Calculate the relative coral cover composed of juveniles over time and functional group.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `is_juvenile` : A boolean vector indicating which size classes are juvenile.
- `location_area` : Vector containing the area of each location with dimensions [locations].

# Returns
Array with dimensions [timesteps ⋅ groups].
"""
function relative_taxa_juveniles(
    relative_cover::AbstractArray{T,4},
    is_juvenile::AbstractVector{Bool},
    location_area::AbstractVector{T}
)::Array{T,2} where {T<:AbstractFloat}
    n_tsteps, n_groups, _, _ = size(relative_cover)
    out_relative_taxa_juveniles::Array{T,2} = zeros(T, n_tsteps, n_groups)
    relative_taxa_juveniles!(
        relative_cover, is_juvenile, location_area, out_relative_taxa_juveniles
    )

    return out_relative_taxa_juveniles
end

"""
    relative_loc_taxa_juveniles!(relative_cover::AbstractArray{T,4}, is_juvenile::AbstractVector{Bool}, location_area::AbstractVector{T}, out_relative_loc_taxa_juveniles::AbstractArray{T,2})::Nothing where {T<:AbstractFloat}

Calculate the relative coral cover copmosed of juveniles over time, functional group and
location. Returns the output into a preallocated buffer

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `is_juvenile` : A boolean vector indicating which size classes are juvenile.
- `out_relative_taxa_juveniles` : Output array buffer with dimensions [timesteps ⋅ groups ⋅ locations].
"""
function relative_loc_taxa_juveniles!(
    relative_cover::AbstractArray{T,4},
    is_juvenile::AbstractVector{Bool},
    out_relative_loc_taxa_juveniles::AbstractArray{T,3}
)::Nothing where {T<:AbstractFloat}
    _is_juveniles = reshape(is_juvenile, (1, 1, :, 1))
    out_relative_loc_taxa_juveniles .=
        dropdims(sum(
                relative_cover .* _is_juveniles; dims=(3,)
    ); dims=(3,))

    return nothing
end

"""
    relative_loc_taxa_juveniles(relative_cover::AbstractArray{T,4}, is_juvenile::AbstractVector{Bool},)::Array{T,3} where {T<:AbstractFloat}

Calculate the relative coral cover composed of juveniles over time, functional group and
location.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `is_juvenile` : A boolean vector indicating which size classes are juvenile.

# Returns
Array with dimensions [timesteps ⋅ groups ⋅ locations].
"""
function relative_loc_taxa_juveniles(
    relative_cover::AbstractArray{T,4},
    is_juvenile::AbstractVector{Bool}
)::Array{T,3} where {T<:AbstractFloat}
    n_tsteps, n_groups, _, n_locs = size(relative_cover)
    out_relative_loc_taxa_juveniles::Array{T,3} = zeros(T, n_tsteps, n_groups, n_locs)
    relative_loc_taxa_juveniles!(
        relative_cover, is_juvenile, out_relative_loc_taxa_juveniles
    )

    return out_relative_loc_taxa_juveniles
end

"""
    absolute_juveniles!(relative_cover::AbstractArray{T,4}, is_juvenile::AbstractVector{Bool}, location_area::AbstractVector{T}, out_absolute_juveniles::AbstractVector{T})::Nothing where {T<:AbstractFloat}

Calculate the coral cover composed of juvenile corals in absolute units for each timesteps.
Write results into a preallocated buffer.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `is_juvenile` : A boolean vector indicating which size classes are juvenile.
- `location_area` : Vector containing the area of each location with dimensions [locations]
- `out_absolute_juveniles` : Output array buffer with dimensions [timesteps].
"""
function absolute_juveniles!(
    relative_cover::AbstractArray{T,4},
    is_juvenile::AbstractVector{Bool},
    location_area::AbstractVector{T},
    out_absolute_juveniles::AbstractVector{T}
)::Nothing where {T<:AbstractFloat}
    _is_juveniles = reshape(is_juvenile, (1, 1, :, 1))
    _location_area = reshape(location_area, (1, 1, 1, :))
    out_absolute_juveniles .=
        dropdims(sum(
                relative_cover .* _location_area .* _is_juveniles; dims=(2, 3, 4)
            ); dims=(2, 3, 4))

    return nothing
end

"""
    absolute_juveniles(relative_cover::AbstractArray{T,4}, is_juvenile::AbstractVector{Bool}, location_area::AbstractVector{T})::Array{T,2} where {T<:AbstractFloat}

Calculate the coral cover composed of juvenile corals in absolute units for each timestep.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `is_juvenile` : A boolean vector indicating which size classes are juvenile.
- `location_area` : Vector containing the area of each location with dimensions [locations]

# Returns
Array buffer with dimensions [timesteps].
"""
function absolute_juveniles(
    relative_cover::AbstractArray{T,4},
    is_juvenile::AbstractVector{Bool},
    location_area::AbstractVector{T}
)::Array{T} where {T<:AbstractFloat}
    n_tsteps = size(relative_cover, 1)
    out_absolute_juveniles::Vector{T} = zeros(T, n_tsteps)
    absolute_juveniles!(relative_cover, is_juvenile, location_area, out_absolute_juveniles)

    return out_absolute_juveniles
end

"""
    absolute_loc_juveniles!(relative_loc_juveniles::AbstractArray{T,2}, k_area::AbstractVector{T}, out_absolute_loc_juveniles::AbstractArray{T,2})::Nothing where {T<:Real}

Calculate the absolute coral cover composed of juveniles for each timestep and location.

# Arguments
- `relative_loc_juveniles` : Relative juvenile cover with dimensions [timesteps ⋅ locations].
- `k_area` : Habitable area for each location.
- `out_absolute_loc_juveniles` : Output array buffer with dimensions [timesteps ⋅ locations].
"""
function absolute_loc_juveniles!(
    relative_cover::AbstractArray{T,4},
    is_juvenile::AbstractVector{Bool},
    location_area::AbstractVector{T},
    out_absolute_loc_juveniles::AbstractArray{T,2}
)::Nothing where {T<:AbstractFloat}
    relative_juveniles!(relative_cover, is_juvenile, out_absolute_loc_juveniles)
    out_absolute_loc_juveniles .*= location_area'

    return nothing
end

"""
    absolute_loc_juveniles(relative_loc_juveniles::AbstractArray{T,2}, k_area::AbstractVector{T})::AbstractArray{T,2} where {T<:Real}

Calculate the absolute coral cover composed of juveniles for each timestep and location.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `is_juvenile` : Boolean mask indicating juvenile size classes.
- `location_area` : Habitable area for each location.

# Returns
A 2D array of absolute juvenile cover with dimensions [timesteps ⋅ locations].
"""
function absolute_loc_juveniles(
    relative_cover::AbstractArray{T,4},
    is_juvenile::AbstractVector{Bool},
    location_area::AbstractVector{T}
)::AbstractArray{T,2} where {T<:AbstractFloat}
    n_timesteps, _, _, n_locations = size(relative_cover)
    if length(location_area) != n_locations
        throw(
            DimensionMismatch(
                "The number of locations in relative_cover and k_area must match."
            )
        )
    end
    out_absolute_loc_juveniles = zeros(T, n_timesteps, n_locations)
    absolute_loc_juveniles!(relative_cover, is_juvenile, location_area, out_absolute_loc_juveniles)

    return out_absolute_loc_juveniles
end

"""
    absolute_taxa_juveniles!(relative_cover::AbstractArray{T,4}, is_juvenile::AbstractVector{Bool}, location_area::AbstractVector{T}, out_absolute_taxa_juveniles::AbstractArray{T,2})::Nothing

Calculate the coral cover occupied by juveniles over timesteps and functional groups. Write
results in a preallocated buffer.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `is_juvenile` : Boolean mask indicating juvenile size classes.
- `location_area` : Habitable area for each location.
- `out_absolute_taxa_juveniles` : Out array buffer with dimensions [timesteps ⋅ groups ⋅ locations]
"""
function absolute_taxa_juveniles!(
    relative_cover::AbstractArray{T,4},
    is_juvenile::AbstractVector{Bool},
    location_area::AbstractVector{T},
    out_absolute_taxa_juveniles::AbstractArray{T,2}
)::Nothing where {T<:AbstractFloat}
    _is_juveniles = reshape(is_juvenile, (1, 1, :, 1))
    _location_area = reshape(location_area, (1, 1, 1, :))
    out_absolute_taxa_juveniles .=
        dropdims(sum(
                relative_cover .* _location_area .* _is_juveniles; dims=(3, 4)
            ); dims=(3, 4))

    return nothing
end

"""
    absolute_taxa_juveniles(relative_cover::AbstractArray{T,4}, is_juvenile::AbstractVector{Bool}, location_area::AbstractVector{T})::Array{T,2} where {T<:AbstractFloat}

Calculate the coral cover occupied by juveniles in absolute units over timesteps, functional
groups and locations.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `is_juvenile` : Boolean mask indicating juvenile size classes.
- `location_area` : Habitable area for each location.

# Returns
A 2D array with dimensions [timesteps ⋅ groups]
"""
function absolute_taxa_juveniles(
    relative_cover::AbstractArray{T,4},
    is_juvenile::AbstractVector{Bool},
    location_area::AbstractVector{T}
)::Array{T,2} where {T<:AbstractFloat}
    n_locs, n_groups, _, _ = size(relative_cover)
    out_absolute_taxa_juveniles::Array{T,2} = zeros(T, n_locs, n_groups)
    absolute_taxa_juveniles!(
        relative_cover, is_juvenile, location_area, out_absolute_taxa_juveniles
    )

    return out_absolute_taxa_juveniles
end

"""
    absolute_loc_taxa_juveniles!(relative_cover::AbstractArray{T,4}, is_juvenile::AbstractVector{Bool}, location_area::AbstractVector{T}, out_absolute_loc_taxa_juveniles::AbstractArray{T,3})::Nothing

Calculate the coral cover occupied by juveniles in absolute units over timesteps, functional
groups and locations. Write results into a preallocated buffer.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `is_juvenile` : Boolean mask indicating juvenile size classes.
- `location_area` : Habitable area for each location.
- `out_absolute_loc_taxa_juveniles` : Out array buffer with dimensions [timesteps ⋅ groups ⋅ locations]
"""
function absolute_loc_taxa_juveniles!(
    relative_cover::AbstractArray{T,4},
    is_juvenile::AbstractVector{Bool},
    location_area::AbstractVector{T},
    out_absolute_loc_taxa_juveniles::AbstractArray{T,3}
)::Nothing where {T<:AbstractFloat}
    _is_juveniles = reshape(is_juvenile, (1, 1, :, 1))
    _location_area = reshape(location_area, (1, 1, 1, :))
    out_absolute_loc_taxa_juveniles .=
        dropdims(sum(
                relative_cover .* _location_area .* _is_juveniles; dims=(3,)
            ); dims=(3,))

    return nothing
end

"""
    absolute_loc_taxa_juveniles(relative_cover::AbstractArray{T,4}, is_juvenile::AbstractVector{Bool}, location_area::AbstractVector{T})::Array{T,3} where {T<:AbstractFloat}

Calculate the coral cover occupied by juveniles in absolute units over timesteps, functional
groups and locations.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `is_juvenile` : Boolean mask indicating juvenile size classes.
- `location_area` : Habitable area for each location.

# Returns
A 2D array of absolute juvenile cover with dimensions [timesteps ⋅ locations].
"""
function absolute_loc_taxa_juveniles(
    relative_cover::AbstractArray{T,4},
    is_juvenile::AbstractVector{Bool},
    location_area::AbstractVector{T}
)::Array{T,3} where {T<:AbstractFloat}
    n_tsteps, n_groups, _, n_locs = size(relative_cover)
    out_absolute_loc_taxa_juveniles = zeros(T, n_tsteps, n_groups, n_locs)
    absolute_loc_taxa_juveniles!(
        relative_cover, is_juvenile, location_area, out_absolute_loc_taxa_juveniles
    )

    return out_absolute_loc_taxa_juveniles
end

"""
    _max_juvenile_area(max_juv_colony_area::T, max_juv_density::T)::T where {T<:Real}

Calculate the maximum possible area that can be covered by juveniles for a given m².

# Arguments
- `max_juv_colony_area` : Maximum colony area of a juvenile in m².
- `max_juv_density` : Maximum juvenile density in individuals/m².
"""
function _max_juvenile_area(max_juv_colony_area::T, max_juv_density::T)::T where {T<:Real}
    return max_juv_density * max_juv_colony_area
end

"""
    _maximum_colony_area(size_spec::AbstractArray{T,2})::T where {T<:AbstractFloat}

Calculate the largest colony size given a range of size classes.
"""
function _maximum_colony_area(mean_colony_diams::AbstractArray{T,2})::T where {T<:AbstractFloat}
    max_idx = argmax(mean_colony_diams)
    mean_colony_size::T = (π / 4) * mean_colony_diams[max_idx]

    return mean_colony_size
end

"""
    juvenile_indicator!(relative_cover::AbstractArray{T,4}, is_juvenile::AbstractVector{Bool}, location_area::AbstractVector{T}, max_juv_colony_area::T, max_juv_density::T)::Nothing

Indicator for juvenile density (0 - 1) where 1 indicates the maximum theoretical density for
juveniles have been achieved.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations]
- `is_juvenile` : Boolean mask indicating juvenile size classes.
- `habitable_area` : Available area habitable by coral for each location.
- `mean_colony_diameters` : Mean colony diameter for each group and size class with dimensions [groups ⋅ size classes]
- `max_juv_density` : Maximum density juveniles can occur in.
- `out_juvenile_indicator` : Output array buffer for the juvenile indicator metrics with dimensions [timesteps ⋅ locations]
"""
function juvenile_indicator!(
    relative_cover::AbstractArray{T,4},
    is_juvenile::AbstractVector{Bool},
    habitable_area::AbstractVector{T},
    mean_colony_diameters::AbstractArray{T,2},
    max_juv_density::T,
    out_juvenile_indicator::AbstractArray{T,2}
)::Nothing where {T<:AbstractFloat}
    # Explicit allocation here
    max_col_area::T = _maximum_colony_area(view(mean_colony_diameters, :, is_juvenile))
    abs_juv = absolute_loc_juveniles(relative_cover, is_juvenile, habitable_area)
    max_juv_area::T = _max_juvenile_area(max_col_area, max_juv_density)
    out_juvenile_indicator .= abs_juv ./ (max_juv_area .* habitable_area')

    return nothing
end

"""
    juvenile_indicator(relative_cover::AbstractArray{T,4}, is_juvenile::AbstractVector{Bool}, location_area::AbstractVector{T}, max_juv_colony_area::T, max_juv_density::T)::AbstractArray{T,2} where {T<:AbstractFloat}

Indicator for juvenile density (0 - 1) where 1 indicates the maximum theoretical density for
juveniles have been achieved. The juvenile indicator is defined as follows.

The maximum juvenile colony area ``J_A`` refers to the maximum mean colony area over all juvenile
size classes and all functional groups. Maximum juvenile density ``J_D``refers to the
maximum density juvenile can occur at over all juvenile size classes and functional groups.
Then Juvenile Indicator (I) is given as,

```math
\\begin{align*}
    I(x; J_A, J_D, H_A) = \\frac{A(x; H_A)}{J_A \\cdot J_D \\cdot H_A},
\\end{align*}
```
where ``H_A`` refers to habitable area and ``A(x; H_A)`` refers to absolute juvenile cover.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations]
- `is_juvenile` : Boolean mask indicating which sizes are juveniles.
- `habitable_area` : Available area habitrable by coral for each location.
- `size_spec` : Array containing the size class diameter bounds with dimensions [groups ⋅ size classes + 1]
- `max_juv_density` : Maximum juvenile density for all juvenile size classes and functional groups.

# Returns
A 2D array of juvenile indicators with dimensions [timesteps ⋅ locations]
"""
function juvenile_indicator(
    relative_cover::AbstractArray{T,4},
    is_juvenile::AbstractVector{Bool},
    habitable_area::AbstractVector{T},
    size_spec::AbstractArray{T,2},
    max_juv_density::T
)::AbstractArray{T,2} where {T<:AbstractFloat}
    n_tsteps, _, _, n_locations = size(relative_cover)
    out_juvenile_indicator = zeros(T, n_tsteps, n_locations)
    juvenile_indicator!(
        relative_cover,
        is_juvenile,
        habitable_area,
        size_spec,
        max_juv_density,
        out_juvenile_indicator
    )

    return out_juvenile_indicator
end

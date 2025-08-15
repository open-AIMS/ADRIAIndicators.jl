"""
    _relative_juveniles!(relative_cover::Array{T,4}, is_juvenile::Vector{Bool}, out_relative_juveniles::Array{T,2})::Nothing where {T<:Real}

Calculate the relative juvenile cover.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `is_juvenile` : A boolean vector indicating which size classes are juvenile.
- `out_relative_juveniles` : Output array buffer with dimensions [timesteps ⋅ locations].
"""
function _relative_juveniles!(
    relative_cover::Array{T,4},
    is_juvenile::Vector{Bool},
    out_relative_juveniles::Array{T,2}
)::Nothing where {T<:Real}
    juvenile_cover = relative_cover[:, :, is_juvenile, :]
    out_relative_juveniles .= dropdims(sum(juvenile_cover; dims=(2, 3)); dims=(2, 3))

    return nothing
end

"""
    relative_juveniles(relative_cover::Array{T,4}, is_juvenile::Vector{Bool})::Array{T,2} where {T<:Real}

Calculate the relative juvenile cover.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `is_juvenile` : A boolean vector indicating which size classes are juvenile.

# Returns
A 2D array of relative juvenile cover with dimensions [timesteps ⋅ locations].
"""
function relative_juveniles(
    relative_cover::Array{T,4},
    is_juvenile::Vector{Bool}
)::Array{T,2} where {T<:Real}
    n_timesteps, _, n_sizes, n_locations = size(relative_cover)
    if length(is_juvenile) != n_sizes
        throw(DimensionMismatch("The length of is_juvenile must match the number of size classes in relative_cover."))
    end
    out_relative_juveniles = zeros(T, n_timesteps, n_locations)
    _relative_juveniles!(relative_cover, is_juvenile, out_relative_juveniles)

    return out_relative_juveniles
end

"""
    _absolute_juveniles!(relative_juveniles::Array{T,2}, k_area::Vector{T}, out_absolute_juveniles::Array{T,2})::Nothing where {T<:Real}

Calculate the absolute juvenile cover.

# Arguments
- `relative_juveniles` : Relative juvenile cover with dimensions [timesteps ⋅ locations].
- `k_area` : Habitable area for each location.
- `out_absolute_juveniles` : Output array buffer with dimensions [timesteps ⋅ locations].
"""
function _absolute_juveniles!(
    relative_cover::Array{T,4},
    is_juvenile::Vector{Bool},
    location_area::Vector{T},
    out_absolute_juveniles::Array{T,2}
)::Nothing where {T<:AbstractFloat}
    _relative_juveniles!(relative_cover, is_juvenile, out_absolute_juveniles)
    out_absolute_juveniles .*= location_area'

    return nothing
end

"""
    absolute_juveniles(relative_juveniles::Array{T,2}, k_area::Vector{T})::Array{T,2} where {T<:Real}

Calculate the absolute juvenile cover.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `is_juvenile` : Boolean mask indicating juvenile size classes.
- `location_area` : Habitable area for each location.

# Returns
A 2D array of absolute juvenile cover with dimensions [timesteps ⋅ locations].
"""
function absolute_juveniles(
    relative_cover::Array{T,4},
    is_juvenile::Vector{Bool},
    location_area::Vector{T}
)::Array{T,2} where {T<:AbstractFloat}
    n_timesteps, _, _, n_locations = size(relative_cover)
    if length(location_area) != n_locations
        throw(DimensionMismatch("The number of locations in relative_juveniles and k_area must match."))
    end
    out_absolute_juveniles = zeros(T, n_timesteps, n_locations)
    _absolute_juveniles!(relative_cover, is_juvenile, location_area, out_absolute_juveniles)

    return out_absolute_juveniles
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
    _juvenile_indicator!(relative_cover::Array{T,4}, is_juvenile::Vector{Bool}, location_area::Vector{T}, max_juv_colony_area::T, max_juv_density::T)::Nothing

Indicator for juvenile density (0 - 1) where 1 indicates the maximum theoretical density for
juveniles have been achieved.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations]
- `is_juvenile` : Boolean mask indicating juvenile size classes.
- `habitable_area` : Available area habitable by coral for each location.
- `max_juv_colony_area` : Maximum juvenile colony area for a single juvenile.
- `max_juv_density` : Maximum density juveniles can occur in.
- `out_juvenile_indicator` : Output array buffer for the juvenile indicator metrics with dimensions [timesteps ⋅ locations]
"""
function _juvenile_indicator!(
    relative_cover::Array{T,4},
    is_juvenile::Vector{Bool},
    habitable_area::Vector{T},
    max_juv_colony_area::T,
    max_juv_density::T,
    out_juvenile_indicator::Array{T,2}
)::Nothing where {T<:AbstractFloat}
    # Explicit allocation here
    abs_juv = absolute_juveniles(relative_cover, is_juvenile, habitable_area)
    max_juv_area::T = _max_juvenile_area(max_juv_colony_area, max_juv_density)
    out_juvenile_indicator .= abs_juv ./ (max_juv_area .* habitable_area)

    return nothing
end

"""
    juvenile_indicator(relative_cover::Array{T,4}, is_juvenile::Vector{Bool}, location_area::Vector{T}, max_juv_colony_area::T, max_juv_density::T)::Array{T,2} where {T<:AbstractFloat}

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations]
- `is_juvenile` : Boolean mask indicating which sizes are juveniles.
- `habitable_area` : Available area habitrable by coral for each location.
- `max_juv_colony_area` : Maximum colony area acrosses all juvenile size classes and functional groups.
- `max_juv_density` : Maximum juvenile density for all juvenile size classes and functional groups.

# Returns
A 2D array of juvenile indicators with dimensions [timesteps ⋅ locations]
"""
function juvenile_indicator(
    relative_cover::Array{T,4},
    is_juvenile::Vector{Bool},
    habitable_area::Vector{T},
    max_juv_colony_area::T,
    max_juv_density::T
)::Array{T,2} where {T<:AbstractFloat}
    n_tsteps, _, _, n_locations = size(relative_cover)
    out_juvenile_indicator = zeros(T, n_tsteps, n_locations)
    _juvenile_indicator!(
        relative_cover,
        is_juvenile,
        habitable_area,
        max_juv_colony_area,
        max_juv_density,
        out_juvenile_indicator
    )

    return out_juvenile_indicator
end

"""
    _dimension_mismatch_message(array_name_1::String, array_name_2::String, dims1::Tuple, dims2::Tuple)::String

Construct an informative error message when a discrepency between array dimensions is detected.

# Arguments
- array_name_1 : Name of the first array.
- array_name_2 : Name of the second array.
- dims1 : Shape of the first array.
- dims2 : Shape of the second array.

# Returns
String containing an informative error message.
"""
function _dimension_mismatch_message(
    array_name_1::String,
    array_name_2::String,
    dims1::Tuple,
    dims2::Tuple
)::String
    msg = "\'$(array_name_1)\' and \'$(array_name_2)\' have mismatching dimensions. "
    msg += "\'$(array_name_1)\' and \'$(array_name_2)\' have shapes, $(dims1) and $(dims2) "
    msg += "respectively. Please check the expected shapes and dimensions are correct."

    return msg
end

"""
    _relative_loc_cover!(relative_cover::Array{T,4}, out_relative_cover::Array{T,2})::Nothing where {T<:Real}

Calculate the relative cover per location by summing over groups and size classes.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `out_relative_cover` : Output array buffer for relative cover with dimensions [timesteps ⋅ locations].
"""
function _relative_loc_cover!(
    relative_cover::Array{T,4},
    out_relative_cover::Array{T,2}
)::Nothing where {T<:Real}
    # Sum over groups and sizes
    out_relative_cover .= dropdims(sum(relative_cover; dims=(2, 3)); dims=(2, 3))

    return nothing
end

"""
    relative_cover(relative_cover::Array{T,4})::Array{T,2} where {T<:Real}

Calculate the relative cover by summing over groups and size classes.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].

# Returns
A 2D array of relative cover with dimensions [timesteps ⋅ locations].
"""
function relative_loc_cover(relative_cover::Array{T,4})::Array{T,2} where {T<:Real}
    n_timesteps, _, _, n_locations = size(relative_cover)
    out_relative_cover = zeros(T, n_timesteps, n_locations)
    _relative_loc_cover!(relative_cover, out_relative_cover)

    return out_relative_cover
end

"""
    _relative_taxa_cover!(relative_cover::Array{T,4}, k_area::Vector{T}, out_relative_taxa_cover::Array{T,2})::Nothing where {T<:Real}

Calculate the relative taxa cover, summed up across all locations.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `k_area` : The coral habitable area for each location.
- `out_relative_taxa_cover` : Output array buffer with dimensions [timesteps ⋅ groups].
"""
function _relative_taxa_cover!(
    relative_cover::Array{T,4},
    k_area::Vector{T},
    out_relative_taxa_cover::Array{T,2}
)::Nothing where {T<:Real}
    # Sum over sizes
    group_cover = dropdims(sum(relative_cover; dims=3); dims=3)  # [timesteps, groups, locations]
    absolute_group_cover = group_cover .* k_area'
    total_k_area = sum(k_area)
    out_relative_taxa_cover .= dropdims(sum(absolute_group_cover; dims=3); dims=3) ./ total_k_area

    return nothing
end

"""
    relative_taxa_cover(relative_cover::Array{T,4}, k_area::Vector{T})::Array{T,2} where {T<:Real}

Calculate the relative taxa cover, summed up across all locations.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `k_area` : The coral habitable area for each location.

# Returns
A 2D array of relative taxa cover with dimensions [timesteps ⋅ groups].
"""
function relative_taxa_cover(
    relative_cover::Array{T,4},
    k_area::Vector{T}
)::Array{T,2} where {T<:Real}
    n_timesteps, n_groups, _, n_locations = size(relative_cover)
    if length(k_area) != n_locations
        throw(DimensionMismatch("The number of locations in relative_cover and k_area must match."))
    end
    out_relative_taxa_cover = zeros(T, n_timesteps, n_groups)
    _relative_taxa_cover!(relative_cover, k_area, out_relative_taxa_cover)

    return out_relative_taxa_cover
end

"""
    _relative_cover!(relative_cover::Array{T,4}, location_area::Vector{T}, out_relative_cover::Vector{T})::Nothing where {T<:AbstractFloat}

For each timestep, calculate the proportion of the entire study area covered by coral.

# Arguments
- `relative_cover` : Raw relative coral cover of dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations]
- `location_areas` : Location area for each location.
- `out_relative_cover` : Vector of dimensions [timesteps].
"""
function _relative_cover!(
    relative_cover::Array{T,4},
    location_area::Vector{T},
    out_relative_cover::Vector{T}
)::Nothing where {T<:AbstractFloat}
    total_area::T = sum(location_area)
    out_relative_cover .= dropdims(sum(
        relative_cover .* reshape(location_area, (1, 1, 1, -1)), dims=(2, 3, 4)
    ), dims=(2, 3, 4)) ./ total_area

    return nothing
end

"""
    relative_cover(relative_cover::Array{T,4}, location_area::Vector{T})::Vector{T} where {T<:AbstractFloat}

For each timestep, calculate the proportion of the entire study area covered by coral.

# Arguments
- `relative_cover` : Raw relative coral cover of dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations]
- `location_areas` : Location area for each location.

# Returns
Vector of dimensions [timesteps].
"""
function relative_cover(
    relative_cover::Array{T,4},
    location_area::Vector{T}
)::Vector{T} where {T<:AbstractFloat}
    n_timesteps = size(relative_cover, 4)
    out_relative_cover = zeros(T, n_timesteps)
    _relative_cover!(relative_cover, location_area, out_relative_cover)

    return out_relative_cover
end

"""
    _relative_loc_taxa_cover!(relative_cover::Array{T,4}, out_relative_loc_taxa_cover::Array{T,3})::Nothing where {T<:Real}

Calculate the relative taxa cover for each location.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `out_relative_loc_taxa_cover` : Output array buffer with dimensions [timesteps ⋅ groups ⋅ locations].
"""
function _relative_loc_taxa_cover!(
    relative_cover::Array{T,4},
    out_relative_loc_taxa_cover::Array{T,3}
)::Nothing where {T<:Real}
    out_relative_loc_taxa_cover .= dropdims(sum(relative_cover; dims=3); dims=3)

    return nothing
end

"""
    relative_loc_taxa_cover(relative_cover::Array{T,4})::Array{T,3} where {T<:Real}

Calculate the relative taxa cover for each location.

# Arguments
- `relative_cover` : Relative cover with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].

# Returns
A 3D array of relative taxa cover with dimensions [timesteps ⋅ groups ⋅ locations].
"""
function relative_loc_taxa_cover(
    relative_cover::Array{T,4}
)::Array{T,3} where {T<:Real}
    n_timesteps, n_groups, _, n_locations = size(relative_cover)
    out_relative_loc_taxa_cover = zeros(T, n_timesteps, n_groups, n_locations)
    _relative_loc_taxa_cover!(relative_cover, out_relative_loc_taxa_cover)

    return out_relative_loc_taxa_cover
end

"""
    _coral_diversity(r_taxa_cover::Array{T, 3}, out_coral_diversity::Array{T,2})::Nothing where {T<:Real}
using Base: _dim_stack

Calculates coral taxa diversity as a dimensionless metric.

# Arguments
- `rel_cover` : Relative Taxa Cover of dimensions [timesteps ⋅ groups ⋅ locations]
- `out_coral_diversity` : Output array buffer [timesteps ⋅ locations]
"""
function _coral_diversity!(
    r_taxa_cover::Array{T,3},
    out_coral_diversity::Array{T,2}
)::Nothing where {T<:Real}
    loc_cover = dropdims(sum(r_taxa_cover; dims=2); dims=2)

    for loc in axes(loc_cover, 2)
        out_coral_diversity[:, loc] =
            1 .- sum((r_taxa_cover[:, :, loc] ./ loc_cover[:, loc]) .^ 2; dims=2)
    end

    replace!(
        out_coral_diversity, NaN => 0.0, Inf => 0.0
    )

    return nothing
end

"""
    coral_diversity(rel_cover::Array{T, 3})::Array{T,2} where {T<:Real}

Calculates coral taxa diversity as a dimensionless metric. Derived from the simpsons diversity.

Formulated as part of a reef condition index by Dr Mike Williams (mjmcwilliam@outlook.com) and
Dr Morgan Pratchett (morgan.pratchett@jcu.edu.au).

The coral diversity metric (``D``) for a given location and timestep is given as

```math
\\begin{aligned}
D(x) = 1 - \\sum_{g=1}^{G} (\\frac{x_g}{x_T})^2,
\\end{aligned}
```

where ``x_g`` is the relative coral cover for the functional group, ``g``, and ``x_T`` is
total relative coral cover at the given location and timestep.

# Arguments
- `rel_cover` : Relative Taxa Cover of dimensions [timesteps ⋅ groups ⋅ locations]

# Returns
Matrix containing coral diversity metric of dimension [timesteps ⋅ locations]
"""
function coral_diversity(rel_cover::Array{T,3})::Array{T,2} where {T<:Real}
    n_tsteps, n_groups, n_locs = size(rel_cover)
    coral_div::Array{T,2} = zeros(T, n_tsteps, n_locs)
    _coral_diversity!(rel_cover, coral_div)

    return coral_div
end

"""
    _coral_evenness!(r_taxa_cover::AbstractArray{T,3}, out_coral_evenness::Array{T,2})::Nothing where {T<:Real}

Calculates evenness across functional coral groups in ADRIA as a diversity metric.
Inverse Simpsons diversity indicator.

# Arguments
- `rel_cover` : Relative Taxa Cover of dimensions [timesteps ⋅ groups ⋅ locations]
- `out_coral_evenness` : Output array buffer [timesteps ⋅ locations]

# References
1. Hill, M. O. (1973).
Diversity and Evenness: A Unifying Notation and Its Consequences.
Ecology, 54(2), 427-432.
https://doi.org/10.2307/1934352
"""
function _coral_evenness!(
    rel_cover::AbstractArray{T,3},
    out_coral_evenness::Array{T,2}
)::Nothing where {T<:Real}
    _, n_grps, _ = size(rel_cover)

    # Sum across groups represents functional diversity
    # Group evenness (Hill 1973, Ecology 54:427-432)
    loc_cover = dropdims(sum(rel_cover; dims=2); dims=2)
    for loc in axes(loc_cover, 2)
        out_coral_evenness[:, loc] =
            1.0 ./ sum((rel_cover[:, :, loc] ./ loc_cover[:, loc]) .^ 2; dims=2)
    end

    out_coral_evenness = replace!(
        out_coral_evenness, NaN => 0.0, Inf => 0.0
    ) ./ n_grps

    return nothing
end

"""
    coral_evenness(r_taxa_cover::AbstractArray{T})::AbstractArray{T} where {T<:Real}

Calculates evenness across functional coral groups in ADRIA as a diversity metric.
Inverse Simpsons diversity indicator.

The coral evenness metric (E) is given as follows,

```math
\\begin{align}
E(x) = \\left(\\sum_{g=1}^{G}\\left(\\frac{x_g}{x_T} \\right)^2\\right)^{-1}
\\end{align}
```

# Arguments
- `rel_cover` : Relative Taxa Cover of dimensions [timesteps ⋅ groups ⋅ locations]

# Returns
Matrix containing coral evenness metric of dimensions [timesteps ⋅ locations]

# References
1. Hill, M. O. (1973).
Diversity and Evenness: A Unifying Notation and Its Consequences.
Ecology, 54(2), 427-432.
https://doi.org/10.2307/1934352
"""
function coral_evenness(rel_cover::Array{T,3})::Array{T,2} where {T<:Real}
    n_steps, _, n_locs = size(rel_cover)
    coral_even::Array{T,2} = zeros(T, n_steps, n_locs)
    _coral_evenness!(rel_cover, coral_even)

    return coral_even
end

"""
    _colony_Lcm2_to_m3m2(colony_mean_diams_cm::Array{T,2}, colony_)::Tuple

Helper function to convert coral colony values from Litres/cm² to m³/m²

# Arguments
- `colony_mean_diams_cm` : Matrix of mean colony diameters
- `planar_area_a` : planar area parameter a
- `planar_area_b` : planar area parameter b

# Returns
Tuple : Assumed colony volume (m³/m²) for each species/size class, theoretical maximum for each species class

# References
1. Aston Eoghan A., Duce Stephanie, Hoey Andrew S., Ferrari Renata (2022).
    A Protocol for Extracting Structural Metrics From 3D Reconstructions of Corals.
    Frontiers in Marine Science, 9.
    https://doi.org/10.3389/fmars.2022.854395

"""
function _colony_Lcm2_to_m3m2(
    colony_mean_area_cm::T,
    planar_area_a::T,
    planar_area_b::T
)::T where {T<:AbstractFloat}
    colony_litres_per_cm2::T = exp(
        planar_area_a + planar_area_b * log(colony_mean_area_cm)
    )
    cm2_to_m3_per_m2::T = 10^-3
    colony_vol_m3_per_m2::T = colony_litres_per_cm2 .* cm2_to_m3_per_m2

    return colony_vol_m3_per_m2
end

"""
    _absolute_shelter_volume!(rel_cover::Array{T,3}, colony_mean_area_cm::Array{T,2}, planar_area_params::Array{T,3}, habitable_area::T, ASV::Array{T,3})::Nothing where {T<:AbstractFloat}

# Arguments
- `rel_cover` : 4-D Array of relative coral cover with dimensions [timesteps ⋅ groups ⋅ size ⋅ locations]
- `colony_mean_area_cm` : Matrix of mean colony diameter with dimensions [groups ⋅ size]
- `planar_area_params` : 3-D array of planar area params with dimensions [groups ⋅ size ⋅ param_type]
- `habitable_area_m2` : Vector of habitable area for each location [locations]
- `out_ASV` : Output array buffer for absolute shelter volume [timesteps ⋅ groups ⋅ size ⋅ locations]
"""
function _absolute_shelter_volume!(
    rel_cover::Array{T,4},
    colony_mean_area_cm::Array{T,2},
    planar_area_params::Array{T,3},
    habitable_area::Vector{T},
    out_ASV::Array{T,4}
)::Nothing where {T<:AbstractFloat}
    colony_vol_m3_per_m2 = _colony_Lcm2_to_m3m2.(
        colony_mean_area_cm,
        view(planar_area_params[:, :, 1]),
        view(planar_area_params[:, :, 2])
    )
    n_groups, n_sizes = size(colony_mean_area_cm)

    abs_cover = rel_cover .* reshape(habitable_area(1, 1, 1, -1))
    out_ASV .= abs_cover .* reshape(colony_vol_m3_per_m2, (1, n_groups, n_sizes, 1))

    return nothing
end

"""
    absolute_shelter_volume(rel_cover::Array{T,4}, colony_mean_area_cm::Array{T,2}, planar_area_params::Array{T,3}, habitable_area::Vector{T})::Array{T,4} where {T<:AbstractFloat}

# Arguments
- `rel_cover` : 4-D Array of relative coral cover with dimensions [timesteps ⋅ groups ⋅ size ⋅ locations]
- `colony_mean_area_cm` : Matrix of mean colony diameter with dimensions [groups ⋅ size]
- `planar_area_params` : 3-D array of planar area params with dimensions [groups ⋅ size ⋅ param_type]
- `habitable_area_m2` : Vector of habitable area for each location [locations]

# Returns
- Output array containing absolute shelter volume [timesteps ⋅ groups ⋅ size ⋅ locations]
"""
function absolute_shelter_volume(
    rel_cover::Array{T,4},
    colony_mean_area_cm::Array{T,2},
    planar_area_params::Array{T,3},
    habitable_area::Vector{T}
)::Array{T,4} where {T<:AbstractFloat}
    out_ASV::Array{T,4} = zeros(T, size(rel_cover)...)
    _absolute_shelter_volume!(
        rel_cover, colony_mean_area_cm, planar_area_params, habitable_area, out_ASV
    )

    return out_ASV
end

"""
    _relative_shelter_volume!(rel_cover::Array{T,4}, colony_mean_area_cm::Array{T,2}, planar_area_params::Array{T,3}, habitable_area_m²::Vector{T}, out_RSV::Array{T,4})::Nothing where {T<:AbstractFloat}

Calculate the relative shelter volume for a range of covers. Relative shelter volume (RSV) is
given by

```math
\\begin{align}
    \\text{RSV}(x) = \frac{ASV(x)}{MSV(x)},
\\end{align}
```

where ASV and MSV are Absolute Shelter Volume and Maximum Shelter Volume respectively.

# Arguments
- rel_cover : Relative Cover array with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- colony_mean_area_cm : Mean colony area per group and size class with dimensions [groups ⋅ sizes].
- planar_area_params : Array containing the planar area parameters with dimensions [groups ⋅ sizes ⋅ param_type].
- habitable_area_m² : Habitable area in m² with dimensions [locations].
- out_RSV : Output Relative shelter volume array buffer with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
"""
function _relative_shelter_volume!(
    rel_cover::Array{T,4},
    colony_mean_area_cm::Array{T,2},
    planar_area_params::Array{T,3},
    habitable_area_m²::Vector{T},
    out_RSV::Array{T,4}
)::Nothing where {T<:AbstractFloat}
    colony_vol_m³_per_m² = _colony_Lcm2_to_m3m2.(
        colony_mean_area_cm,
        view(planar_area_params[:, :, 1]),
        view(planar_area_params[:, :, 2])
    )

    n_groups::Int64, n_sizes::Int64 = size(colony_mean_area_cm)
    n_locations::Int64 = length(habitable_area_m²)

    abs_cover_m²::Array{T,4} = rel_cover .* reshape(habitable_area_m², (1, 1, 1, -1))
    ASV_m³ = abs_cover_m² .* reshape(colony_vol_m³_per_m², (1, n_groups, n_sizes, 1))

    max_colony_vol_m³::Vector{T} = dropdims(maximum(colony_vol_m³_per_m², dims=2), dims=2)
    # Calculate maximum shelter volume m³ [group ⋅ location]
    MSV_m³::Matrix{T} = habitable_area_m²' .* reshape(max_colony_vol_m³, (-1, 1))
    out_RSV .= ASV_m³ ./ reshape(MSV_m³, (1, n_groups,1 , n_locations))

    return nothing
end

"""
    relative_shelter_volume(relative_cover::Array{T,4}, colony_mean_area_cm::Array{T,2}, planar_area_params::Array{T,3}, habitable_area_m²::Vector{T})::Array{T,4} where {T<:Real}

Calculate the relative shelter volume for a range of covers. Relative shelter volume (RSV) is
given by

```math
\\begin{align}
    \\text{RSV}(x) = \\frac{ASV(x)}{MSV(x)},
\\end{align}
```

where ASV and MSV are Absolute Shelter Volume and Maximum Shelter Volume respectively.

# Arguments
- rel_cover : Relative Cover array with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- colony_mean_area_cm : Mean colony area per group and size class with dimensions [groups ⋅ sizes].
- planar_area_params : Array containing the planar area parameters with dimensions [groups ⋅ sizes ⋅ param_type].
- habitable_area_m² : Habitable area in m² with dimensions [locations].

# Returns
Relative shelter volume in an array with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations]
"""
function relative_shelter_volume(
    relative_cover::Array{T,4},
    colony_mean_area_cm::Array{T,2},
    planar_area_params::Array{T,3},
    habitable_area_m²::Vector{T}
)::Array{T,4} where {T<:Real}
    n_tsteps::Int64, n_groups::Int64, n_sizes::Int64, n_locs::Int64 = size(relative_cover)

    if size(colony_mean_area_cm) != (n_groups, n_sizes)
        throw(DimensionMismatch(_dimension_mismatch_message(
            "relative_cover",
            "colony_mean_area_cm",
            size(relative_cover),
            size(colony_mean_area_cm)
        )))
    end
    if size(planar_area_params) != (n_groups, n_sizes, 2)
        throw(DimensionMismatch(_dimension_mismatch_message(
            "relative_cover",
            "planar_area_params",
            size(relative_cover),
            size(planar_area_params)
        )))
    end
    if size(habitable_area_m²) != n_locs
        throw(DimensionMismatch(_dimension_mismatch_message(
            "relative_cover",
            "habitable_area_m²",
            size(relative_cover),
            size(habitable_area_m²)
        )))
    end

    RSV::Array{T,4} = zeros(T, n_tsteps, n_groups, n_sizes, n_locs)
    _relative_shelter_volume!(
        relative_cover, colony_mean_area_cm, planar_area_params, habitable_area_m², RSV
    )

    return RSV
end

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
    relative_juveniles::Array{T,2},
    k_area::Vector{T},
    out_absolute_juveniles::Array{T,2}
)::Nothing where {T<:Real}
    out_absolute_juveniles .= relative_juveniles .* k_area'
    return nothing
end

"""
    absolute_juveniles(relative_juveniles::Array{T,2}, k_area::Vector{T})::Array{T,2} where {T<:Real}

Calculate the absolute juvenile cover.

# Arguments
- `relative_juveniles` : Relative juvenile cover with dimensions [timesteps ⋅ locations].
- `k_area` : Habitable area for each location.

# Returns
A 2D array of absolute juvenile cover with dimensions [timesteps ⋅ locations].
"""
function absolute_juveniles(
    relative_juveniles::Array{T,2},
    k_area::Vector{T}
)::Array{T,2} where {T<:Real}
    n_timesteps, n_locations = size(relative_juveniles)
    if length(k_area) != n_locations
        throw(DimensionMismatch("The number of locations in relative_juveniles and k_area must match."))
    end
    out_absolute_juveniles = zeros(T, n_timesteps, n_locations)
    _absolute_juveniles!(relative_juveniles, k_area, out_absolute_juveniles)
    return out_absolute_juveniles
end

"""
    _max_juvenile_area(max_colony_area_m2::T, max_juv_density::T)::T where {T<:Real}

Calculate the maximum possible area that can be covered by juveniles for a given m².

# Arguments
- `max_colony_area_m2` : Maximum colony area of a juvenile in m².
- `max_juv_density` : Maximum juvenile density in individuals/m².
"""
function _max_juvenile_area(max_colony_area_m2::T, max_juv_density::T)::T where {T<:Real}
    return max_juv_density * max_colony_area_m2
end

"""
    _juvenile_indicator!(absolute_juveniles::Array{T,2}, k_area::Vector{T}, max_colony_area_m2::T, max_juv_density::T, out_juvenile_indicator::Array{T,2})::Nothing where {T<:Real}

Calculate the juvenile indicator.

# Arguments
- `absolute_juveniles` : Absolute juvenile cover with dimensions [timesteps ⋅ locations].
- `k_area` : Habitable area for each location.
- `max_colony_area_m2` : Maximum colony area of a juvenile in m².
- `max_juv_density` : Maximum juvenile density in individuals/m².
- `out_juvenile_indicator` : Output array buffer with dimensions [timesteps ⋅ locations].
"""
function _juvenile_indicator!(
    absolute_juveniles::Array{T,2},
    k_area::Vector{T},
    max_colony_area_m2::T,
    max_juv_density::T,
    out_juvenile_indicator::Array{T,2}
)::Nothing where {T<:Real}
    max_juv_area = _max_juvenile_area(max_colony_area_m2, max_juv_density)
    usable_k_area = max.(k_area, 1.0)'
    out_juvenile_indicator .= absolute_juveniles ./ (max_juv_area .* usable_k_area)
    return nothing
end

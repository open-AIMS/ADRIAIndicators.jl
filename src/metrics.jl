"""
    _dimension_mismatch_message(array_name_1::String, array_name_2::String, dims1::Tuple, dims2::Tuple)::String

Construct an informative error message when a discrepancy between array dimensions is detected.

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
    _coral_diversity(r_taxa_cover::Array{T, 3}, out_coral_diversity::Array{T,2})::Nothing where {T<:Real}

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

Calculates coral taxa diversity as a dimensionless metric. Derived from the Simpson's Diversity Index.

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
- `planar_area_int` : planar area model intercept
- `planar_area_coef` : planar area model coefficient

# Returns
Tuple : Assumed colony volume (m³/m²) for each species/size class

# References
1. Aston Eoghan A., Duce Stephanie, Hoey Andrew S., Ferrari Renata (2022).
    A Protocol for Extracting Structural Metrics From 3D Reconstructions of Corals.
    Frontiers in Marine Science, 9.
    https://doi.org/10.3389/fmars.2022.854395

"""
function _colony_Lcm2_to_m3m2(
    colony_mean_area_cm::T,
    planar_area_int::T,
    planar_area_coef::T
)::T where {T<:AbstractFloat}
    colony_litres_per_cm2::T = exp(
        planar_area_int + planar_area_coef * log(colony_mean_area_cm)
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
- `planar_area_params` : 3-D array of planar area parameters with dimensions [groups ⋅ size ⋅ (intercept, coefficient)]
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
        view(planar_area_params,:,:,1),
        view(planar_area_params,:,:,2)
    )
    n_groups, n_sizes = size(colony_mean_area_cm)

    abs_cover = rel_cover .* reshape(habitable_area, (1, 1, 1, :))
    out_ASV .= abs_cover .* reshape(colony_vol_m3_per_m2, (1, n_groups, n_sizes, 1))

    return nothing
end

"""
    absolute_shelter_volume(rel_cover::Array{T,4}, colony_mean_area_cm::Array{T,2}, planar_area_params::Array{T,3}, habitable_area::Vector{T})::Array{T,4} where {T<:AbstractFloat}

Calculate the volume of shelter provided by the given coral cover. This function uses
log-log linear models to predict the volume of shelter provided from a given planar area.
The parametrisation of this log-log linear model must be provided by the user. The log-log
linear model is given by

```math
\\begin{align}
    \\log(S) = b + a\\log(PA),
\\end{align}
```
where ``S`` and ``PA`` are shelter volume (``m^3m_{-2}``) and planar area, respectively.
Then absolute shelter volume is given by

```math
\\begin{align}
    ASV = A_C \\cdot S,
\\end{align}
where ``ASV`` and ``A_C`` refers to absolute shelter volume and absolute coral cover,
respectively.
```

# Arguments
- `rel_cover` : 4-D Array of relative coral cover with dimensions [timesteps ⋅ groups ⋅ size ⋅ locations]
- `colony_mean_area_cm` : Matrix of mean colony diameter with dimensions [groups ⋅ size]
- `planar_area_params` : 3-D array of planar area params with dimensions [groups ⋅ size ⋅ (intercept, coefficient)]
- `habitable_area_m2` : Vector of habitable area for each location [locations]

# Returns
- Output array containing absolute shelter volume [timesteps ⋅ groups ⋅ size ⋅ locations]

# References
1. Urbina-Barreto, I., Chiroleu, F., Pinel, R., Fréchon, L., Mahamadaly, V.,
     Elise, S., Kulbicki, M., Quod, J.-P., Dutrieux, E., Garnier, R.,
     Henrich Bruggemann, J., Penin, L., & Adjeroud, M. (2021).
   Quantifying the shelter capacity of coral reefs using photogrammetric
     3D modeling: From colonies to reefscapes.
   Ecological Indicators, 121, 107151.
   https://doi.org/10.1016/j.ecolind.2020.107151
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

Calculate the relative shelter volume for a range of covers. Relative to the theoretical
maximum of 50% cover of a coral species with the largest colony volume.
Relative shelter volume (RSV) is given by

```math
\\begin{align}
    \\text{RSV}(x) = \frac{ASV(x)}{MSV(x)},
\\end{align}
```

where ASV and MSV are Absolute Shelter Volume and Maximum Shelter Volume respectively.

# Arguments
- rel_cover : Relative Cover array with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- colony_mean_area_cm : Mean colony area per group and size class with dimensions [groups ⋅ sizes].
- planar_area_params : Array containing the planar area parameters with dimensions [groups ⋅ sizes ⋅ (intercept, coefficient)].
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

    abs_cover_m²::Array{T,4} = rel_cover .* reshape(habitable_area_m², (1, 1, 1, :))
    ASV_m³ = abs_cover_m² .* reshape(colony_vol_m³_per_m², (1, n_groups, n_sizes, 1))

    max_colony_vol_m³::T = max(colony_vol_m³_per_m²)
    # Calculate maximum shelter volume m³ [group ⋅ location]
    MSV_m³::Vector{T} = habitable_area_m²' .* max_colony_vol_m³ .* 0.5
    out_RSV .= ASV_m³ ./ reshape(MSV_m³, (1, 1, 1, n_locations))

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
- planar_area_params : Array containing the planar area parameters with dimensions [groups ⋅ sizes ⋅ (intercept, coefficient)].
- habitable_area_m² : Habitable area in m² with dimensions [locations].

# Returns
Relative shelter volume in an array with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations]

# References
1. Urbina-Barreto, I., Chiroleu, F., Pinel, R., Fréchon, L., Mahamadaly, V.,
     Elise, S., Kulbicki, M., Quod, J.-P., Dutrieux, E., Garnier, R.,
     Henrich Bruggemann, J., Penin, L., & Adjeroud, M. (2021).
   Quantifying the shelter capacity of coral reefs using photogrammetric
     3D modeling: From colonies to reefscapes.
   Ecological Indicators, 121, 107151.
   https://doi.org/10.1016/j.ecolind.2020.107151
"""
function relative_shelter_volume(
    relative_cover::Array{T,4},
    colony_mean_area_cm::Array{T,2},
    planar_area_params::Array{T,3},
    habitable_area_m²::Vector{T}
)::Array{T,4} where {T<:Real}
    n_tsteps::Int64, n_groups::Int64, n_sizes::Int64, n_locs::Int64 = size(relative_cover)

    if size(colony_mean_area_cm) != (n_groups, n_sizes)
        throw(
            DimensionMismatch(
                _dimension_mismatch_message(
                    "relative_cover",
                    "colony_mean_area_cm",
                    size(relative_cover),
                    size(colony_mean_area_cm)
                )
            )
        )
    end
    if size(planar_area_params) != (n_groups, n_sizes, 2)
        throw(
            DimensionMismatch(
                _dimension_mismatch_message(
                    "relative_cover",
                    "planar_area_params",
                    size(relative_cover),
                    size(planar_area_params)
                )
            )
        )
    end
    if size(habitable_area_m²) != n_locs
        throw(
            DimensionMismatch(
                _dimension_mismatch_message(
                    "relative_cover",
                    "habitable_area_m²",
                    size(relative_cover),
                    size(habitable_area_m²)
                )
            )
        )
    end

    RSV::Array{T,4} = zeros(T, n_tsteps, n_groups, n_sizes, n_locs)
    _relative_shelter_volume!(
        relative_cover, colony_mean_area_cm, planar_area_params, habitable_area_m², RSV
    )

    return RSV
end

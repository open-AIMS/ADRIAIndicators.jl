"""
    _dimension_mismatch_message(array_name_1::String, array_name_2::String, dims1::Tuple, dims2::Tuple)::String

Construct an informative error message when a discrepancy between array dimensions is detected.

# Arguments
- `array_name_1` : Name of the first array.
- `array_name_2` : Name of the second array.
- `dims1` : Shape of the first array.
- `dims2` : Shape of the second array.

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
    msg *= "\'$(array_name_1)\' and \'$(array_name_2)\' have shapes, $(dims1) and $(dims2) "
    msg *= "respectively. Please check the expected shapes and dimensions are correct."

    return msg
end

"""
    coral_diversity(r_taxa_cover::Array{T, 3}, out_coral_diversity::Array{T,2})::Nothing where {T<:Real}

Calculates coral taxa diversity as a dimensionless metric.

# Arguments
- `r_taxa_cover` : Relative Taxa Cover of dimensions [timesteps ⋅ groups ⋅ locations], relative to habitable area.
- `out_coral_diversity` : Output array buffer [timesteps ⋅ locations]
"""
function coral_diversity!(
    r_taxa_cover::AbstractArray{T,3},
    out_coral_diversity::Array{T,2}
)::Nothing where {T<:Real}
    n_tsteps, n_groups, n_locs = size(r_taxa_cover)

    for l in 1:n_locs
        for t in 1:n_tsteps
            loc_cover = zero(T)
            for g in 1:n_groups
                loc_cover += r_taxa_cover[t, g, l]
            end

            if loc_cover > 0.0
                sum_sq = zero(T)
                for g in 1:n_groups
                    sum_sq += (r_taxa_cover[t, g, l] / loc_cover)^2
                end
                out_coral_diversity[t, l] = 1.0 - sum_sq
            else
                out_coral_diversity[t, l] = 0.0
            end
        end
    end

    return nothing
end

"""
    coral_diversity(rel_cover::Array{T, 3})::Array{T,2} where {T<:Real}

Calculates coral taxa diversity as a dimensionless metric. Derived from the Simpson's Diversity Index.

Formulated as part of a reef condition index by Dr Mike Williams (mjmcwilliam@outlook.com) and
Dr Morgan Pratchett (morgan.pratchett@jcu.edu.au).

The coral diversity metric (``D``) for a given location and timestep is given as

```math
\\begin{align*}
D(x) = 1 - \\sum_{g=1}^{G} (\\frac{x_g}{x_T})^2,
\\end{align*}
```

where ``x_g`` is the relative coral cover for the functional group, ``g``, and ``x_T`` is
total relative coral cover at the given location and timestep.

# Arguments
- `rel_cover` : Relative Taxa Cover of dimensions [timesteps ⋅ groups ⋅ locations], relative to habitable area.

# Returns
Matrix containing coral diversity metric of dimension [timesteps ⋅ locations]
"""
function coral_diversity(rel_cover::Array{T,3})::Array{T,2} where {T<:Real}
    n_tsteps, _, n_locs = size(rel_cover)
    coral_div::Array{T,2} = zeros(T, n_tsteps, n_locs)
    coral_diversity!(rel_cover, coral_div)

    return coral_div
end

"""
    coral_evenness!(r_taxa_cover::AbstractArray{T,3}, out_coral_evenness::Array{T,2})::Nothing where {T<:Real}

Calculates evenness across functional coral groups in ADRIA as a diversity metric.
Inverse Simpsons diversity indicator.

# Arguments
- `rel_cover` : Relative Taxa Cover of dimensions [timesteps ⋅ groups ⋅ locations], relative to habitable area.
- `out_coral_evenness` : Output array buffer [timesteps ⋅ locations]

# References
1. Hill, M. O. (1973).
    Diversity and Evenness: A Unifying Notation and Its Consequences.
    Ecology, 54(2), 427-432.
    https://doi.org/10.2307/1934352
"""
function coral_evenness!(
    rel_cover::AbstractArray{T,3},
    out_coral_evenness::Array{T,2}
)::Nothing where {T<:Real}
    n_tsteps, n_groups, n_locs = size(rel_cover)

    for l in 1:n_locs
        for t in 1:n_tsteps
            loc_cover = zero(T)
            for g in 1:n_groups
                loc_cover += rel_cover[t, g, l]
            end

            if loc_cover > 0.0
                sum_sq = zero(T)
                for g in 1:n_groups
                    sum_sq += (rel_cover[t, g, l] / loc_cover)^2
                end
                out_coral_evenness[t, l] = 1.0 / sum_sq
            else
                out_coral_evenness[t, l] = 0.0
            end
        end
    end

    return nothing
end

"""
    coral_evenness(r_taxa_cover::AbstractArray{T})::AbstractArray{T} where {T<:Real}

Calculates evenness across functional coral groups in ADRIA as a diversity metric.
Inverse Simpsons diversity indicator.

The coral evenness metric (E) is given as follows,

```math
\\begin{align*}
E(x) = \\left(\\sum_{g=1}^{G}\\left(\\frac{x_g}{x_T} \\right)^2\\right)^{-1}
\\end{align*}
```

# Arguments
- `rel_cover` : Relative Taxa Cover of dimensions [timesteps ⋅ groups ⋅ locations], relative to habitable area.

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
    coral_evenness!(rel_cover, coral_even)

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
    colony_mean_diam_cm::T,
    planar_area_int::T,
    planar_area_coef::T
)::T where {T<:AbstractFloat}
    colony_litres_per_cm2::T = exp(
        planar_area_int + planar_area_coef * log(colony_mean_diam_cm)
    )
    dm3_to_m3_per_m2::T = 10^-3
    colony_vol_m3_per_m2::T = colony_litres_per_cm2 .* dm3_to_m3_per_m2

    return colony_vol_m3_per_m2
end

"""
    maximum_colony_volume(colony_mean_diam_cm::AbstractArray{T,2}, planar_area_params::AbstractArray{T,3})::T where {T<:AbstractFloat}

Find the maximum colony volume per m² among all provided species and size classes.
"""
function maximum_colony_volume(
    colony_mean_diam_cm::AbstractArray{T,2},
    planar_area_params::AbstractArray{T,3}
)::T where {T<:AbstractFloat}
    max_vol::T = zero(T)
    for idx in CartesianIndices(colony_mean_diam_cm)
        v = _colony_Lcm2_to_m3m2(
            colony_mean_diam_cm[idx],
            planar_area_params[idx, 1],
            planar_area_params[idx, 2]
        )
        if v > max_vol
            max_vol = v
        end
    end

    return max_vol
end

"""
    maximum_colony_volume(diam::T, intercept::T, coefficient::T)::T where {T<:AbstractFloat}

Calculate colony volume per m² for a specific parameterization.
"""
function maximum_colony_volume(
    diam::T,
    intercept::T,
    coefficient::T
)::T where {T<:AbstractFloat}
    return _colony_Lcm2_to_m3m2(diam, intercept, coefficient)
end

"""
    absolute_shelter_volume!(rel_cover::AbstractArray{T,3}, colony_mean_diam_cm::AbstractArray{T,2}, planar_area_params::AbstractArray{T,3}, habitable_area::T, ASV::AbstractArray{T,3})::Nothing where {T<:AbstractFloat}

# Arguments
- `rel_cover` : 4-D Array of relative coral cover with dimensions [timesteps ⋅ groups ⋅ size ⋅ locations], relative to habitable area.
- `colony_mean_diam_cm` : Matrix of mean colony diameter with dimensions [groups ⋅ size]
- `planar_area_params` : 3-D array of planar area parameters with dimensions [groups ⋅ size ⋅ (intercept, coefficient)]
- `habitable_area_m2` : Vector of habitable area for each location [locations]
- `out_ASV` : Output array buffer for absolute shelter volume [timesteps ⋅ groups ⋅ size ⋅ locations]
"""
function absolute_shelter_volume!(
    rel_cover::AbstractArray{T,4},
    colony_mean_diam_cm::AbstractArray{T,2},
    planar_area_params::AbstractArray{T,3},
    habitable_area::AbstractVector{T},
    out_ASV::AbstractArray{T,4}
)::Nothing where {T<:AbstractFloat}
    n_groups, n_sizes = size(colony_mean_diam_cm)
    n_timesteps, _, _, n_locations = size(rel_cover)

    for l in 1:n_locations
        h_area = habitable_area[l]
        for s in 1:n_sizes
            for g in 1:n_groups
                colony_vol = _colony_Lcm2_to_m3m2(
                    colony_mean_diam_cm[g, s],
                    planar_area_params[g, s, 1],
                    planar_area_params[g, s, 2]
                )
                for t in 1:n_timesteps
                    out_ASV[t, g, s, l] = rel_cover[t, g, s, l] * h_area * colony_vol
                end
            end
        end
    end

    return nothing
end

"""
    absolute_shelter_volume(rel_cover::AbstractArray{T,4}, colony_mean_diam_cm::AbstractArray{T,2}, planar_area_params::AbstractArray{T,3}, habitable_area::AbstractVector{T})::AbstractArray{T,4} where {T<:AbstractFloat}

Calculate the volume of shelter provided by the given coral cover. This function uses
log-log linear models to predict the volume of shelter provided from a given planar area.
The parametrisation of this log-log linear model must be provided by the user. For
possible parametrisations of the log-log linear model used to
predict shelter volume from planar area, see Urbina-Barreto et al. (2021), [1].

The log-log
linear model is given by

```math
\\begin{align*}
    \\log(S) = b + a\\log(PA),
\\end{align*}
```
where ``S`` and ``PA`` are shelter volume (``m^3m_{-2}``) and planar area, respectively.
Then absolute shelter volume is given by

```math
\\begin{align*}
    ASV = A_C \\cdot S,
\\end{align*}
```
where ``ASV`` and ``A_C`` refers to absolute shelter volume and absolute coral cover,
respectively.

# Arguments
- `rel_cover` : 4-D Array of relative coral cover with dimensions [timesteps ⋅ groups ⋅ size ⋅ locations], relative to habitable area.
- `colony_mean_diam_cm` : Matrix of mean colony diameter with dimensions [groups ⋅ size]
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
    rel_cover::AbstractArray{T,4},
    colony_mean_diam_cm::AbstractArray{T,2},
    planar_area_params::AbstractArray{T,3},
    habitable_area::AbstractVector{T}
)::Array{T,4} where {T<:AbstractFloat}
    out_ASV::Array{T,4} = zeros(T, size(rel_cover)...)
    absolute_shelter_volume!(
        rel_cover, colony_mean_diam_cm, planar_area_params, habitable_area, out_ASV
    )

    return out_ASV
end

"""
    relative_shelter_volume!(rel_cover::AbstractArray{T,4}, colony_mean_diam_cm::AbstractArray{T,2}, planar_area_params::AbstractArray{T,3}, habitable_area_m²::AbstractVector{T}, out_RSV::AbstractArray{T,4}, reference::Tuple{T, T, T})::Nothing where {T<:AbstractFloat}

Calculate the relative shelter volume for a range of covers. Relative to the theoretical
maximum of 50% cover of a coral species with the specified reference parameterisation.
Relative shelter volume (RSV) is given by

```math
\\begin{align*}
    \\text{RSV}(x) = \\frac{ASV(x)}{MSV(x)},
\\end{align*}
```

where ASV and MSV are Absolute Shelter Volume and Maximum Shelter Volume respectively.

# Arguments
- `rel_cover` : Relative Cover array with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations], relative to habitable area.
- `colony_mean_diam_cm` : Mean colony diameter per group and size class with dimensions [groups ⋅ sizes].
- `planar_area_params` : Array containing the planar area parameters with dimensions [groups ⋅ sizes ⋅ (intercept, coefficient)].
- `habitable_area_m²` : Habitable area in m² with dimensions [locations].
- `out_RSV` : Output Relative shelter volume array buffer with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations].
- `reference` : Parameterisation to use as reference (mean diameter, intercept, coefficient).
"""
function relative_shelter_volume!(
    rel_cover::AbstractArray{T,4},
    colony_mean_diam_cm::AbstractArray{T,2},
    planar_area_params::AbstractArray{T,3},
    habitable_area_m²::AbstractVector{T},
    out_RSV::AbstractArray{T,4},
    reference::Tuple{T, T, T}
)::Nothing where {T<:AbstractFloat}
    n_groups, n_sizes = size(colony_mean_diam_cm)
    n_timesteps, _, _, n_locations = size(rel_cover)

    # Calculate max colony volume per m² from reference
    max_colony_vol::T = _colony_Lcm2_to_m3m2(reference...)

    for l in 1:n_locations
        # Maximum shelter volume m³ = habitable_area * max_colony_vol * 0.5
        # Since we want relative shelter volume:
        # RSV = (rel_cover * habitable_area * colony_vol) / (habitable_area * max_colony_vol * 0.5)
        # RSV = (rel_cover * colony_vol) / (max_colony_vol * 0.5)
        msv_factor = max_colony_vol * 0.5
        for s in 1:n_sizes
            for g in 1:n_groups
                colony_vol = _colony_Lcm2_to_m3m2(
                    colony_mean_diam_cm[g, s],
                    planar_area_params[g, s, 1],
                    planar_area_params[g, s, 2]
                )
                factor = colony_vol / msv_factor
                for t in 1:n_timesteps
                    out_RSV[t, g, s, l] = rel_cover[t, g, s, l] * factor
                end
            end
        end
    end

    return nothing
end

"""
    relative_shelter_volume(relative_cover::AbstractArray{T,4}, colony_mean_diam_cm::AbstractArray{T,2}, planar_area_params::AbstractArray{T,3}, habitable_area_m²::AbstractVector{T}, reference::Tuple{T, T, T})::AbstractArray{T,4} where {T<:Real}

Calculate the relative shelter volume for a range of covers. Relative shelter volume (RSV) is
given by

```math
\\begin{align*}
    \\text{RSV}(x) = \\frac{ASV(x)}{MSV(x)},
\\end{align*}
```

where ASV and MSV are Absolute Shelter Volume and Maximum Shelter Volume respectively.
The maximum shelter volume is defined by assuming the maximum theoretical shelter volume
produced by a specified reference parameterisation at 50% cover.
For possible parametrisations of the log-log linear model used to predict shelter volume
from planar area, see Urbina-Barreto et al., [1].

# Arguments
- `rel_cover` : Relative Cover array with dimensions [timesteps ⋅ groups ⋅ sizes ⋅ locations], relative to habitable area.
- `colony_mean_diam_cm` : Mean colony area per group and size class with dimensions [groups ⋅ sizes].
- `planar_area_params` : Array containing the planar area parameters with dimensions [groups ⋅ sizes ⋅ (intercept, coefficient)].
- `habitable_area_m²` : Habitable area in m² with dimensions [locations].
- `reference` : Parameterisation to use as reference (mean diameter, intercept, coefficient).

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
    relative_cover::AbstractArray{T,4},
    colony_mean_diam_cm::AbstractArray{T,2},
    planar_area_params::AbstractArray{T,3},
    habitable_area_m²::AbstractVector{T},
    reference::Tuple{T, T, T}
)::Array{T,4} where {T<:Real}
    n_tsteps::Int64, n_groups::Int64, n_sizes::Int64, n_locs::Int64 = size(relative_cover)

    if size(colony_mean_diam_cm) != (n_groups, n_sizes)
        throw(
            DimensionMismatch(
                _dimension_mismatch_message(
                    "relative_cover",
                    "colony_mean_diam_cm",
                    size(relative_cover),
                    size(colony_mean_diam_cm)
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
    if size(habitable_area_m², 1) != n_locs
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
    relative_shelter_volume!(
        relative_cover,
        colony_mean_diam_cm,
        planar_area_params,
        habitable_area_m²,
        RSV,
        reference
    )

    return RSV
end

"""
    scenario_metric(metric::AbstractArray{T}, location_area::AbstractVector{T}, location_dim::Int; is_relative::Bool=true, return_relative::Bool=true)::AbstractArray{T} where {T<:AbstractFloat}

Aggregate a metric across the location dimension.

This function can take a metric that is either relative to location area or absolute, and
can return a metric that is either relative to the total area or absolute.

# Arguments
- `metric` : An array containing the metric to aggregate.
- `location_area` : A vector of area values for each location.
- `location_dim` : The dimension of the `metric` array that corresponds to location.
- `is_relative` : Whether the input `metric` is relative to location area. Defaults to `true`.
- `return_relative` : Whether the output should be relative to total area. Defaults to `true`.

# Returns
An array with the location dimension removed, containing the aggregated metric.
"""
function scenario_metric(
    metric::AbstractArray{T}, location_area::AbstractVector{T}, location_dim::Int;
    is_relative::Bool=true, return_relative::Bool=true
)::AbstractArray{T} where {T<:AbstractFloat}
    dims = ones(Int, ndims(metric))
    dims[location_dim] = length(location_area)
    _location_area = reshape(location_area, dims...)

    absolute_metric = is_relative ? metric .* _location_area : metric
    aggregated_metric = sum(absolute_metric; dims=location_dim)

    total_area = sum(location_area)
    result = return_relative ? aggregated_metric ./ total_area : aggregated_metric

    return dropdims(result; dims=location_dim)
end

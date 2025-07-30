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

Calculates coral taxa diversity as a dimensionless metric. Derived from the simpsons diversity.

Formulated as part of a reef condition index by Dr Mike Williams (mjmcwilliam@outlook.com) and
Dr Morgan Pratchett (morgan.pratchett@jcu.edu.au).

The coral diversity metric (``D``) for a given location and timestep if given as

```math
\\begin{aligned}
D(x) = \\sum_{g=1}^{G} (\\frac{x_g}{x_T})^2,
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
    coral_div::Array{T,2} = zeros(T, n_steps, n_locs)
    _coral_diversity!(rel_cover, coral_div)

    return coral_div
end

"""
    coral_evenness(r_taxa_cover::AbstractArray{T})::AbstractArray{T} where {T<:Real}

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

    replace!(
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
E(x) = \\left(\sum_{g=1}^{G}\\left(\\frac{x_g}{x_T} \\right)^2\\right)^{-1}
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

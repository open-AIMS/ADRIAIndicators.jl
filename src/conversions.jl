const AF = AbstractFloat

"""
    relative_cover_to_ltmp_cover!(relative_cover::Array{T,N}, habitable_area_m²::AbstractVector{T}, reef_area_m²::AbstractVector{T}, location_dim::Int64, out_ltmp_cover::Array{T,N})::Nothing where {T<:AbstractFloat,N}

Convert relative cover to LTMP cover. If location dimension is -1, then assume the location
dimension has already been aggregated and multiple by the ratio of total reef area to total
habitable area.

# Arguments
- `relative_cover` : Relative cover (relative to habitable area) with value in [0, 1].
- `habitable_area_m²` : Habitable area of reef in m².
- `reef_area_m²` : Area of entire reef in m².
- `location_dim` : Index of the location dimensions. Use -1 for no location dimension.
- `out_ltmp_cover` : Array buffer of the same shape as relative_cover.
"""
function relative_cover_to_ltmp_cover!(
    relative_cover::AbstractArray{<:AF,N},
    habitable_area_m²::AbstractVector{<:AF},
    reef_area_m²::AbstractVector{<:AF},
    location_dim::Int64,
    out_ltmp_cover::Array{<:AF,N}
)::Nothing where {N}
    if location_dim == -1
        out_ltmp_cover .= relative_cover .* sum(habitable_area_m²) ./ sum(reef_area_m²)
        return nothing
    end
    dims = ones(Int, ndims(relative_cover))
    dims[location_dim] = length(habitable_area_m²)
    reshape_dims = Tuple(dims)

    area_coefficient = reshape(habitable_area_m² ./ reef_area_m², reshape_dims)
    out_ltmp_cover .= relative_cover .* area_coefficient

    return nothing
end

"""
    relative_cover_to_ltmp_cover(relative_cover::Array{T}, habitable_area_m²::AbstractVector{T}, reef_area_m²::AbstractVector{T}, location_dim::Int64)::Array{T} where {T<:AbstractFloat}

Convert relative cover to LTMP cover. The conversion is given by

```math
\\begin{align*}
\\text{LTMP} = \\text{RC} \\cdot \\frac{A_H}{A_R},
\\end{align*}
```

where LTMP, RC, ``A_H`` and ``A_R`` represent LTMP cover, relative cover,
habitable area and reef area respectively. 

If location dimension is -1, then assume the location
dimension has already been aggregated and multiple by the ratio of total reef area to total
habitable area.

# Arguments
- `relative_cover` : Relative cover (relative to habitable area) with value in [0, 1].
- `habitable_area_m²` : Habitable area of reef in m².
- `reef_area_m²` : Area of entire reef in m².
- `location_dim` : Index of the location dimensions. Use -1 for no location dimension.

# Returns
LTMP cover with same array shape as the input `relative_cover`.
"""
function relative_cover_to_ltmp_cover(
    relative_cover::AbstractArray{<:AF},
    habitable_area_m²::AbstractVector{<:AF},
    reef_area_m²::AbstractVector{<:AF},
    location_dim::Int64
)::Array{<:AF} 
    out_ltmp_cover = zeros(eltype(relative_cover), size(relative_cover)...)
    relative_cover_to_ltmp_cover!(
        relative_cover,
        habitable_area_m²,
        reef_area_m²,
        location_dim,
        out_ltmp_cover
    )

    return out_ltmp_cover
end

"""
    ltmp_cover_to_relative_cover!(ltmp_cover::Array{T,N}, habitable_area_m²::AbstractVector{T}, reef_area_m²::AbstractVector{T}, location_dim::Int64, out_relative_cover::Array{T,N})::Nothing where {T<:AbstractFloat,N}

Convert LTMP cover to relative cover. If location dimension is -1, then assume the location
dimension has already been aggregated and multiple by the ratio of total reef area to total
habitable area.

# Arguments
- `ltmp_cover` : LTMP cover with value in [0, 1].
- `habitable_area_m²` : Habitable area of reef in m².
- `reef_area_m²` : Area of entire reef in m².
- `location_dim` : Index of the location dimensions. Use -1 for no location dimension.
- `out_relative_cover` : Array buffer of the same shape as ltmp_cover.
"""
function ltmp_cover_to_relative_cover!(
    ltmp_cover::AbstractArray{<:AF,N},
    habitable_area_m²::AbstractVector{<:AF},
    reef_area_m²::AbstractVector{<:AF},
    location_dim::Int64,
    out_relative_cover::Array{<:AF,N}
)::Nothing where {N}
    if location_dim == -1
        out_relative_cover .= relative_cover .* sum(reef_area_m²) ./ sum(habitable_area_m²)
        return nothing
    end
    dims = ones(Int, ndims(ltmp_cover))
    dims[location_dim] = length(habitable_area_m²)
    reshape_dims = Tuple(dims)

    area_coefficient = reshape(reef_area_m² ./ habitable_area_m², reshape_dims)
    out_relative_cover .= ltmp_cover .* area_coefficient

    return nothing
end

"""
    ltmp_cover_to_relative_cover(ltmp_cover::Array{T}, habitable_area_m²::AbstractVector{T}, reef_area_m²::AbstractVector{T}, location_dim::Int64)::Array{T} where {T<:AbstractFloat}

Convert LTMP cover to relative cover. The conversion is given by

```math
\\begin{align*}
\\text{RC} = \\text{LTMP} \\cdot \\frac{A_R}{A_H},
\\end{align*}
```

where LTMP, RC, ``A_H`` and ``A_R`` represent LTMP cover, relative cover,
habitable area and reef area respectively. 

If location dimension is -1, then assume the location
dimension has already been aggregated and multiple by the ratio of total reef area to total
habitable area.


# Arguments
- `ltmp_cover` : LTMP cover with value in [0, 1].
- `habitable_area_m²` : Habitable area of reef in m².
- `reef_area_m²` : Area of entire reef in m².
- `location_dim` : Index of the location dimensions. Use -1 for no location dimension.

# Returns
Relative cover with same array shape as the input `ltmp_cover`.
"""
function ltmp_cover_to_relative_cover(
    ltmp_cover::AbstractArray{<:AF},
    habitable_area_m²::AbstractVector{<:AF},
    reef_area_m²::AbstractVector{<:AF},
    location_dim::Int64
)::Array{<:AF}
    out_relative_cover = zeros(eltype(ltmp_cover), size(ltmp_cover)...)
    ltmp_cover_to_relative_cover!(
        ltmp_cover,
        habitable_area_m²,
        reef_area_m²,
        location_dim,
        out_relative_cover
    )

    return out_relative_cover
end

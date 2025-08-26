"""
    _rhc_to_rrc!(relative_habitable_cover::Array{T}, habitable_area_m²::AbstractVector{T}, reef_area_m²::AbstractVector{T}, location_dim::Int64, out_rrc::Array{T,N})::Nothing where {T<:AbstractFloat}

Convert relative habitable cover to relative reef cover.

# Arguments
- `relative_reef_cover` : Relative reef cover with value in [0, 1]
- `habitable_area_m²` : Habitable area of reef.
- `reef_area_m²` : Area of entire reef.
- `location_dim` : Index of the location dimensions. For example location_dim=3 if the 
third dimension in `relative_reef_cover` is the location dimensions.
- `out_rrc` : Array buffer of the same shape as relative_reef_cover
"""
function rhc_to_rrc!(
    relative_habitable_cover::Array{T,N},
    habitable_area_m²::AbstractVector{T},
    reef_area_m²::AbstractVector{T},
    location_dim::Int64,
    out_rrc::Array{T,N}
)::Nothing where {T<:AbstractFloat,N}
    reshape_idx::Tuple = Tuple(
        i == location_dim ? -1 : 1 for i in 1:ndims(relative_habitable_cover)
    )
    area_coefficient::Array{T} = reshape(habitable_area_m² ./ reef_area_m², reshape_idx)
    out_rrc .= relative_habitable_cover .* area_coefficient

    return nothing
end

"""
    rhc_to_rrc(relative_habitable_cover::Array{T}, habitable_area_m²::AbstractVector{T}, reef_area_m²::AbstractVector{T}, location_dim::Int64)::Array{T} where {T<:AbstractFloat}

Convert relative habitable cover to relative reef cover. The conversion is given by

```math
\\begin{align*}
\\text{RRC} = \\text{RHC} \\cdot \\frac{A_H}{A_R},
\\end{align*}
```

where RRC, RHC, ``A_H`` and ``A_R`` represent relative reef cover, relative habitable cover
habitable area and reef area respectively.

# Arguments
- `relative_reef_cover` : Relative reef cover with value in [0, 1]
- `habitable_area_m²` : Habitable area of reef.
- `reef_area_m²` : Area of entire reef.
- `location_dim` : Index of the location dimensions. For example location_dim=3 if the 
third dimension in `relative_reef_cover` is the location dimsnions.

# Returns
Relative reef cover with same array shape as the input `relative_reef_cover`.
"""
function rhc_to_rrc(
    relative_habitable_cover::Array{T},
    habitable_area_m²::AbstractVector{T},
    reef_area_m²::AbstractVector{T},
    location_dim::Int64
)::Array{T} where {T<:AbstractFloat}
    out_rrc::Array{T} = zeros(T, size(relative_habitable_cover)...)
    rhc_to_rrc!(
        relative_habitable_cover,
        habitable_area_m²,
        reef_area_m²,
        location_dim,
        out_rrc
    )

    return out_rrc
end

"""
    rrc_to_rhc!(relative_reef_cover::Array{T}, habitable_area_m²::AbstractVector{T}, reef_area_m²::AbstractVector{T}, location_dim::Int64, out_rhc::Array{T,N})::Nothing where {T<:AbstractFloat}

Convert relative reef cover to relative habitable cover.

# Arguments
- `relative_reef_cover` : Relative reef cover with value in [0, 1]
- `habitable_area_m²` : Habitable area of reef.
- `reef_area_m²` : Area of entire reef.
- `location_dim` : Index of the location dimensions. For example location_dim=3 if the 
third dimension in `relative_reef_cover` is the location dimsnions.
- `out_rhc` : Array buffer of the same shape as relative_reef_cover
"""
function rrc_to_rhc!(
    relative_reef_cover::Array{T,N},
    habitable_area_m²::AbstractVector{T},
    reef_area_m²::AbstractVector{T},
    location_dim::Int64,
    out_rhc::Array{T,N}
)::Array{T} where {T<:AbstractFloat,N}
    reshape_idx::Tuple = Tuple(
        i == location_dim ? -1 : 1 for i in 1:ndims(relative_reef_cover)
    )
    area_coefficient::AbstractVector{T} = reshape(reef_area_m² ./ habitable_area_m², reshape_idx)
    out_rhc .= relative_reef_cover .* area_coefficient

    return nothing
end

"""
    rrc_to_rhc(relative_reef_cover::Array{T}, habitable_area_m²::AbstractVector{T}, reef_area_m²::AbstractVector{T}, location_dim::Int64)::Array{T} where {T<:AbstractFloat}

Convert relative reef cover to relative habitable cover. The conversion is given by

```math
\\begin{align*}
\\text{RHC} = \\text{RRC} \\cdot \\frac{A_R}{A_H},
\\end{align*}
```

where RRC, RHC, ``A_H`` and ``A_R`` represenht relative reef cover, relative habitable cover
habitable area and reef area respectively.


# Arguments
- `relative_reef_cover` : Relative reef cover with value in [0, 1]
- `habitable_area_m²` : Habitable area of reef.
- `reef_area_m²` : Area of entire reef.
- `location_dim` : Index of the location dimensions. For example location_dim=3 if the 
third dimension in `relative_reef_cover` is the location dimsnions.

# Returns
Relative habitable cover with same array shape as the input `relative_reef_cover`.
"""
function rrc_to_rhc(
    relative_reef_cover::Array{T},
    habitable_area_m²::AbstractVector{T},
    reef_area_m²::AbstractVector{T},
    location_dim::Int64
)::Array{T} where {T<:AbstractFloat}
    out_rhc::Array{T,N} = zeros(T, size(relative_reef_cover)...)
    rrc_to_rhc!(
        relative_reef_cover,
        habitable_area_m²,
        reef_area_m²,
        location_dim,
        out_rhc
    )

    return out_rhc
end

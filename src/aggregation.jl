"""
    _aggregate_relative_cover!(relative_cover::Array{T,2}, location_area::Vector{T}, out_aggregate_cover::Vector{T})::Nothing where {T<:AbstractFloat}
"""
function _aggregate_relative_cover!(
    relative_cover::Array{T,2},
    location_area::Vector{T},
    out_aggregate_cover::Vector{T}
)::Nothing where {T<:AbstractFloat}
    total_area::T = sum(location_area)
    out_aggregate_cover .= dropdims(sum(
        relative_cover .* location_area', dims=2
    ), dims=2) ./ total_area

    return nothing
end

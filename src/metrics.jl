"""
    coral_diversity(r_taxa_cover::Array{T, 3})::AbstractArray{T} where {T<:Real}

Calculates coral taxa diversity as a dimensionless metric. Derived from the simpsons diversity,
D = 1-sum_i((cov_i/cov)^2) where cov is the total coral cover and cov_i is the cover of taxa i.
Formulated as part of a reef condition index by Dr Mike Williams (mjmcwilliam@outlook.com) and
Dr Morgan Pratchett (morgan.pratchett@jcu.edu.au).

# Arguments
- r_taxa_cover : Relative Taxa Cover of dimensions [timesteps ⋅ groups ⋅ locations]
"""
function coral_diversity!(
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

Base.@ccallable function coral_diversity(
    n_tsteps::Cint,
    n_groups::Cint,
    n_locs::Cint,
    relative_taxa_cover::Ptr{Float64},
    output_taxa_cover::Ptr{Float64}
)::Cvoid
    wrapped_r_taxa_cover::Array{Float64, 3} = unsafe_wrap(
        Array,
        relative_taxa_cover,
        (n_tsteps, n_groups, n_locs)
    )
    wrapped_output::Array{Float64, 2} = unsafe_wrap(
        Array,
        output_taxa_cover,
        (n_tsteps, n_locs)
    )

    coral_diversity!(wrapped_r_taxa_cover, wrapped_output)
end

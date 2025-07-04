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

"""
    coral_diversity(n_tsteps::Cint, n_groups::Cint, n_locs::Cint, relative_taxa_cover::Ptr{Float64}, output_taxa_cover::Ptr{Float64})

DOCSTRING

# Arguments:
- `n_tsteps`: DESCRIPTION
- `n_groups`: DESCRIPTION
- `n_locs`: DESCRIPTION
- `relative_taxa_cover`: DESCRIPTION
- `output_taxa_cover`: DESCRIPTION
"""
Base.@ccallable function coral_diversity(
    n_tsteps::Cint,
    n_groups::Cint,
    n_locs::Cint,
    relative_taxa_cover::Ptr{Float64},
    output_taxa_cover::Array{Float64, 2}
)::Cvoid
    println("here")
    wrapped_r_taxa_cover::Array{Float64, 3} = unsafe_wrap(
        Array,
        relative_taxa_cover,
        (n_tsteps, n_groups, n_locs)
    )
    #wrapped_output::Array{Float64, 2} = unsafe_wrap(
    #    Array,
    #    output_taxa_cover,
    #    (n_tsteps, n_locs)
    #)

    coral_diversity!(wrapped_r_taxa_cover, output_taxa_cover)#wrapped_output)
end

# """
#     coral_evenness(n_timesteps::Cint, n_groups::Cint, n_locs::Cint, relative_taxa_cover::Ptr{Float64}, output_evenness::Ptr{Float64})
# 
# DOCSTRING
# 
# # Arguments:
# - `n_timesteps`: Number of timesteps in input data arrays.
# - `n_groups`: DESCRIPTION
# - `n_locs`: DESCRIPTION
# - `relative_taxa_cover`: DESCRIPTION
# - `output_evenness`: DESCRIPTION
# """
# Base.@ccallable function coral_evenness(
#     n_timesteps::Cint,
#     n_groups::Cint,
#     n_locs::Cint,
#     relative_taxa_cover::Ptr{Float64},  # Input:  Array [timesteps ⋅ groups ⋅ locations]
#     output_evenness::Ptr{Float64}       # Output: Array [timesteps ⋅ locations]
# )::Cvoid
#     return ""
# end
# 
# """
#     juvenile_indicator(n_timesteps::Cint, n_species::Cint, n_locs::Cint, proportional_cover::Ptr{Float64}, class_ids::Ptr{Cint}, mean_colony_diams::Ptr{Float64}, k_area::Ptr{Float64}, output_juvenile_indicator::Ptr{Float64})
# 
# DOCSTRING
# 
# # Arguments:
# - `n_timesteps`: Number of timesteps in input data arrays.
# - `n_species`: DESCRIPTION
# - `n_locs`: DESCRIPTION
# - `proportional_cover`: DESCRIPTION
# - `class_ids`: DESCRIPTION
# - `mean_colony_diams`: DESCRIPTION
# - `k_area`: DESCRIPTION
# - `output_juvenile_indicator`: DESCRIPTION
# """
# Base.@ccallable function juvenile_indicator(
#     n_timesteps::Cint,
#     n_species::Cint,
#     n_locs::Cint,
#     proportional_cover::Ptr{Float64},   # Input: Array [timesteps ⋅ species ⋅ locations]
#     class_ids::Ptr{Cint},               # Input: Vector [species]
#     mean_colony_diams::Ptr{Float64},    # Input: Vector [species]
#     k_area::Ptr{Float64},               # Input: Vector [locations]
#     output_juvenile_indicator::Ptr{Float64} # Output: Array [timesteps ⋅ locations]
# )::Cvoid
#     return nothing
# end
# 
# """
#     absolute_shelter_volume(n_timesteps::Cint, n_species::Cint, n_locs::Cint, proportional_cover::Ptr{Float64}, k_area::Ptr{Float64}, colony_vol_m3_per_m2::Ptr{Float64}, output_asv::Ptr{Float64})
# 
# DOCSTRING
# 
# # Arguments:
# - `n_timesteps`: Number of timesteps in input data arrays.
# - `n_species`: DESCRIPTION
# - `n_locs`: DESCRIPTION
# - `proportional_cover`: DESCRIPTION
# - `k_area`: DESCRIPTION
# - `colony_vol_m3_per_m2`: DESCRIPTION
# - `output_asv`: DESCRIPTION
# """
# Base.@ccallable function absolute_shelter_volume(
#     n_timesteps::Cint,
#     n_species::Cint,
#     n_locs::Cint,
#     proportional_cover::Ptr{Float64},  # Input: Array [timesteps ⋅ species ⋅ locations]
#     k_area::Ptr{Float64},              # Input: Vector [locations]
#     colony_vol_m3_per_m2::Ptr{Float64},# Input: Vector [species]
#     output_asv::Ptr{Float64}           # Output: Array [timesteps ⋅ locations]
# )::Cvoid
#     return nothing
# end
# 
# """
#     relative_shelter_volume(n_timesteps::Cint, n_species::Cint, n_locs::Cint, n_groups::Cint, proportional_cover::Ptr{Float64}, k_area::Ptr{Float64}, colony_vol_m3_per_m2::Ptr{Float64}, max_colony_vol_m3_per_m2::Ptr{Float64}, output_rsv::Ptr{Float64})
# 
# DOCSTRING
# 
# # Arguments:
# - `n_timesteps`: Number of timesteps in input data arrays.
# - `n_species`: DESCRIPTION
# - `n_locs`: DESCRIPTION
# - `n_groups`: DESCRIPTION
# - `proportional_cover`: DESCRIPTION
# - `k_area`: DESCRIPTION
# - `colony_vol_m3_per_m2`: DESCRIPTION
# - `max_colony_vol_m3_per_m2`: DESCRIPTION
# - `output_rsv`: DESCRIPTION
# """
# Base.@ccallable function relative_shelter_volume(
#     n_timesteps::Cint,
#     n_species::Cint,
#     n_locs::Cint,
#     n_groups::Cint,
#     proportional_cover::Ptr{Float64},    # Input: Array [timesteps ⋅ species ⋅ locations]
#     k_area::Ptr{Float64},                # Input: Vector [locations]
#     colony_vol_m3_per_m2::Ptr{Float64},  # Input: Vector [species]
#     max_colony_vol_m3_per_m2::Ptr{Float64},# Input: Vector [groups]
#     output_rsv::Ptr{Float64}             # Output: Array [timesteps ⋅ locations]
# )::Cvoid
#     return nothing
# end
# 
# """
#     reef_condition_index(n_timesteps::Cint, n_locs::Cint, relative_cover::Ptr{Float64}, evenness::Ptr{Float64}, relative_shelter_volume::Ptr{Float64}, juvenile_indicator::Ptr{Float64}, threshold::Cint, output_rci::Ptr{Float64})
# 
# DOCSTRING
# 
# # Arguments:
# - `n_timesteps`: Number of timesteps in input data arrays.
# - `n_locs`: DESCRIPTION
# - `relative_cover`: DESCRIPTION
# - `evenness`: DESCRIPTION
# - `relative_shelter_volume`: DESCRIPTION
# - `juvenile_indicator`: DESCRIPTION
# - `threshold`: DESCRIPTION
# - `output_rci`: DESCRIPTION
# """
# Base.@ccallable function reef_condition_index(
#     n_timesteps::Cint,
#     n_locs::Cint,
#     relative_cover::Ptr{Float64},           # Input: Array [timesteps ⋅ locations]
#     evenness::Ptr{Float64},                 # Input: Array [timesteps ⋅ locations]
#     relative_shelter_volume::Ptr{Float64},  # Input: Array [timesteps ⋅ locations]
#     juvenile_indicator::Ptr{Float64},       # Input: Array [timesteps ⋅ locations]
#     threshold::Cint,                        # Input: Scalar
#     output_rci::Ptr{Float64}                # Output: Array [timesteps ⋅ locations]
# )::Cvoid
#     return nothing
# end
# 
# """
#     reef_tourism_index(n_timesteps::Cint, n_locs::Cint, relative_cover::Ptr{Float64}, evenness::Ptr{Float64}, relative_shelter_volume::Ptr{Float64}, juvenile_indicator::Ptr{Float64}, intercept_uncertainty::Float64, output_rti::Ptr{Float64})
# 
# DOCSTRING
# 
# # Arguments:
# - `n_timesteps`: Number of timesteps in input data arrays.
# - `n_locs`: DESCRIPTION
# - `relative_cover`: DESCRIPTION
# - `evenness`: DESCRIPTION
# - `relative_shelter_volume`: DESCRIPTION
# - `juvenile_indicator`: DESCRIPTION
# - `intercept_uncertainty`: DESCRIPTION
# - `output_rti`: DESCRIPTION
# """
# Base.@ccallable function reef_tourism_index(
#     n_timesteps::Cint,
#     n_locs::Cint,
#     relative_cover::Ptr{Float64},           # Input: Array [timesteps ⋅ locations]
#     evenness::Ptr{Float64},                 # Input: Array [timesteps ⋅ locations]
#     relative_shelter_volume::Ptr{Float64},  # Input: Array [timesteps ⋅ locations]
#     juvenile_indicator::Ptr{Float64},       # Input: Array [timesteps ⋅ locations]
#     intercept_uncertainty::Float64,         # Input: Scalar
#     output_rti::Ptr{Float64}                # Output: Array [timesteps ⋅ locations]
# )::Cvoid
#     return nothing
# end
# 
# """
#     reef_fish_index(n_timesteps::Cint, n_locs::Cint, relative_cover::Ptr{Float64}, intercept_uncertainty_1::Float64, intercept_uncertainty_2::Float64, output_rfi::Ptr{Float64})
# 
# DOCSTRING
# 
# # Arguments:
# - `n_timesteps`: Number of timesteps in input data arrays.
# - `n_locs`: DESCRIPTION
# - `relative_cover`: DESCRIPTION
# - `intercept_uncertainty_1`: DESCRIPTION
# - `intercept_uncertainty_2`: DESCRIPTION
# - `output_rfi`: DESCRIPTION
# """
# Base.@ccallable function reef_fish_index(
#     n_timesteps::Cint,
#     n_locs::Cint,
#     relative_cover::Ptr{Float64},         # Input: Array [timesteps ⋅ locations]
#     intercept_uncertainty_1::Float64,     # Input: Scalar
#     intercept_uncertainty_2::Float64,     # Input: Scalar
#     output_rfi::Ptr{Float64}              # Output: Array [timesteps ⋅ locations]
# )::Cvoid
#     return nothing
# end
# 
# # ------------------ Non-indices
# 
# """
#     relative_cover(n_timesteps::Cint, n_groups::Cint, n_locs::Cint, relative_cover::Ptr{Float64}, output_relative_cover::Ptr{Float64})
# 
# DOCSTRING
# 
# # Arguments:
# - `n_timesteps`: Number of timesteps in input data arrays.
# - `n_groups`: DESCRIPTION
# - `n_locs`: DESCRIPTION
# - `relative_cover`: DESCRIPTION
# - `output_relative_cover`: DESCRIPTION
# """
# Base.@ccallable function relative_cover(
#     n_timesteps::Cint,
#     n_groups::Cint,
#     n_locs::Cint,
#     relative_cover::Ptr{Float64},     # Input:  Array [timesteps ⋅ groups ⋅ locations]
#     output_relative_cover::Ptr{Float64}   # Output: Array [timesteps ⋅ locations]
# )::Cvoid
#     return nothing
# end
# 
# """
#     total_absolute_cover(n_timesteps::Cint, n_locs::Cint, relative_cover::Ptr{Float64}, k_area::Ptr{Float64}, output_absolute_cover::Ptr{Float64})
# 
# DOCSTRING
# 
# # Arguments:
# - `n_timesteps`: Number of timesteps in input data arrays.
# - `n_locs`: DESCRIPTION
# - `relative_cover`: DESCRIPTION
# - `k_area`: DESCRIPTION
# - `output_absolute_cover`: DESCRIPTION
# """
# Base.@ccallable function total_absolute_cover(
#     n_timesteps::Cint,
#     n_locs::Cint,
#     relative_cover::Ptr{Float64},         # Input: Array [timesteps ⋅ locations]
#     k_area::Ptr{Float64},                 # Input: Vector [locations]
#     output_absolute_cover::Ptr{Float64}   # Output: Array [timesteps ⋅ locations]
# )::Cvoid
#     return nothing
# end
# 
# """
#     relative_taxa_cover(n_timesteps::Cint, n_groups::Cint, n_size_classes::Cint, n_locs::Cint, proportional_cover::Ptr{Float64}, k_area::Ptr{Float64}, output_taxa_cover::Ptr{Float64})
# 
# DOCSTRING
# 
# # Arguments:
# - `n_timesteps`: Number of timesteps in input data arrays.
# - `n_groups`: DESCRIPTION
# - `n_size_classes`: DESCRIPTION
# - `n_locs`: DESCRIPTION
# - `proportional_cover`: DESCRIPTION
# - `k_area`: DESCRIPTION
# - `output_taxa_cover`: DESCRIPTION
# """
# Base.@ccallable function relative_taxa_cover(
#     n_timesteps::Cint,
#     n_groups::Cint,
#     n_size_classes::Cint,
#     n_locs::Cint,
#     proportional_cover::Ptr{Float64},     # Input:  Array [timesteps ⋅ groups ⋅ sizes ⋅ locations]
#     k_area::Ptr{Float64},                 # Input:  Vector [locations]
#     output_taxa_cover::Ptr{Float64}       # Output: Array [timesteps ⋅ groups]
# )::Cvoid
#     return nothing
# end
# 
# """
#     relative_loc_taxa_cover(n_timesteps::Cint, n_species::Cint, n_locs::Cint, n_groups::Cint, proportional_cover::Ptr{Float64}, k_area::Ptr{Float64}, output_loc_taxa_cover::Ptr{Float64})
# 
# DOCSTRING
# 
# # Arguments:
# - `n_timesteps`: Number of timesteps in input data arrays.
# - `n_species`: DESCRIPTION
# - `n_locs`: DESCRIPTION
# - `n_groups`: DESCRIPTION
# - `proportional_cover`: DESCRIPTION
# - `k_area`: DESCRIPTION
# - `output_loc_taxa_cover`: DESCRIPTION
# """
# Base.@ccallable function relative_loc_taxa_cover(
#     n_timesteps::Cint,
#     n_species::Cint,
#     n_locs::Cint,
#     n_groups::Cint,
#     proportional_cover::Ptr{Float64},      # Input:  Array [timesteps ⋅ species ⋅ locations]
#     k_area::Ptr{Float64},                  # Input:  Vector [locations]
#     output_loc_taxa_cover::Ptr{Float64}    # Output: Array [timesteps ⋅ groups ⋅ locations]
# )::Cvoid
#     return nothing
# end
# 
# """
#     relative_juveniles(n_timesteps::Cint, n_groups::Cint, n_locs::Cint, proportional_cover::Ptr{Float64}, class_ids::Ptr{Cint}, output_relative_juveniles::Ptr{Float64})
# 
# DOCSTRING
# 
# # Arguments:
# - `n_timesteps`: Number of timesteps in input data arrays.
# - `n_groups`: DESCRIPTION
# - `n_locs`: DESCRIPTION
# - `proportional_cover`: DESCRIPTION
# - `class_ids`: DESCRIPTION
# - `output_relative_juveniles`: DESCRIPTION
# """
# Base.@ccallable function relative_juveniles(
#     n_timesteps::Cint,
#     n_groups::Cint,
#     n_locs::Cint,
#     proportional_cover::Ptr{Float64},       # Input:  Array [timesteps ⋅ group0s ⋅ locations]
#     class_ids::Ptr{Cint},                   # Input:  Vector [groups]
#     output_relative_juveniles::Ptr{Float64} # Output: Array [timesteps ⋅ locations]
# )::Cvoid
#     return nothing
# end
# 
# """
#     absolute_juveniles(n_timesteps::Cint, n_species::Cint, n_locs::Cint, proportional_cover::Ptr{Float64}, class_ids::Ptr{Cint}, k_area::Ptr{Float64}, output_absolute_juveniles::Ptr{Float64})
# 
# DOCSTRING
# 
# # Arguments:
# - `n_timesteps`: Number of timesteps in input data arrays.
# - `n_species`: DESCRIPTION
# - `n_locs`: DESCRIPTION
# - `proportional_cover`: DESCRIPTION
# - `class_ids`: DESCRIPTION
# - `k_area`: DESCRIPTION
# - `output_absolute_juveniles`: DESCRIPTION
# """
# Base.@ccallable function absolute_juveniles(
#     n_timesteps::Cint,
#     n_species::Cint,
#     n_locs::Cint,
#     proportional_cover::Ptr{Float64},      # Input:  Array [timesteps ⋅ species ⋅ locations]
#     class_ids::Ptr{Cint},                  # Input:  Vector [species]
#     k_area::Ptr{Float64},                  # Input:  Vector [locations]
#     output_absolute_juveniles::Ptr{Float64}# Output: Array [timesteps ⋅ locations]
# )::Cvoid
#     return nothing
# end

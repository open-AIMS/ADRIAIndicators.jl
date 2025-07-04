using ReefMetrics
using Random

n_locs::Int32 = 10
n_tsteps::Int32 = 9
n_groups::Int32 = 5

dummy_input::Array{Float64, 3} = rand(n_tsteps, n_groups, n_locs)
dummy_output::Array{Float64, 2} = zeros(Float64, n_tsteps, n_locs)

ReefMetrics.coral_diversity(
    n_tsteps, 
    n_groups, 
    n_locs, 
    pointer(dummy_input), 
    dummy_output#pointer(dummy_output)
)

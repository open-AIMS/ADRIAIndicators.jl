"""
Tests for shelter volume metrics.
"""

using Test

using ReefMetrics: absolute_shelter_volume, relative_shelter_volume

@testset "Absolute Shelter Volume" begin
    n_tsteps, n_groups, n_sizes, n_locs = 2, 2, 3, 2
    relative_cover = zeros(Float64, n_tsteps, n_groups, n_sizes, n_locs)
    planar_area_params = zeros(Float64, n_groups, n_sizes, 2)
    habitable_area = [100.0, 200.0]

    relative_cover[1, 1, :, 1] = [0.1, 0.1, 0.05]
    relative_cover[1, 2, :, 1] = [0.05, 0.05, 0.2]

    relative_cover[2, 1, :, 2] = [0.0, 0.05, 0.3]
    relative_cover[2, 2, :, 2] = [0.2, 0.0, 0.05]

    colony_mean_area_cm = [
        15.0 20.0 25.0;
        20.0 25.0 30.0;
    ]
    planar_area_params[1, :, :] = [
        -8.32 2.5;
        -7.37 2.34;
        -8.31 2.47
    ]
    planar_area_params[2, :, :] = [
        -7.37 2.34;
        -8.31 2.47;
        -9.69 2.49
    ]

    sv = exp.(planar_area_params[:, :, 1] .+ planar_area_params[:, :, 2] .* log.(colony_mean_area_cm))
    sv .*= 0.001

    asv = absolute_shelter_volume(
        relative_cover, colony_mean_area_cm, planar_area_params, habitable_area
    )
    @test size(asv) == (n_tsteps, n_groups, n_sizes, n_locs)
    @test all(asv .≈ (
        relative_cover .* reshape(habitable_area, (1, 1, 1, :)) .* reshape(sv, (1, 2, 3, 1))
    ))
end

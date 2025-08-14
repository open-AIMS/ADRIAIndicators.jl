"""
This file tests

- all relative cover metrics
- juvenile metrics
- shelter volume metrics
"""

using Test
using ReefMetrics

@testset "Relative Juveniles" begin
    n_tsteps, n_groups, n_sizes, n_locs = 2, 2, 4, 2
    relative_cover = zeros(Float64, n_tsteps, n_groups, n_sizes, n_locs)

    is_juvenile = [true, true, false, false]

    relative_cover[1, 1, :, 1] = [0.1, 0.1, 0.05, 0.05]
    relative_cover[1, 2, :, 1] = [0.05, 0.05, 0.2, 0.2]

    relative_cover[2, 1, :, 2] = [0.0, 0.05, 0.3, 0.3]
    relative_cover[2, 2, :, 2] = [0.2, 0.0, 0.05, 0.0]

    @testset "Normal Cases" begin
        # Test the allocating version
        rel_juv = relative_juveniles(relative_cover, is_juvenile)
        @test size(rel_juv) == (n_tsteps, n_locs)

        @test rel_juv[1, 1] ≈ 0.3
        @test rel_juv[1, 2] == 0.0
        @test rel_juv[2, 1] == 0.0
        @test rel_juv[2, 2] == 0.25
    end

    @testset "Edge Cases" begin
        # Test with no juvenile classes
        no_juveniles = [false, false, false, false]
        @test all(relative_juveniles(relative_cover, no_juveniles) .== 0.0)

        # Test with all classes as juveniles
        all_juveniles = [true, true, true, true]
        expected_sum = dropdims(sum(relative_cover, dims=(2,3)), dims=(2,3))
        @test relative_juveniles(relative_cover, all_juveniles) ≈ expected_sum
    end
end

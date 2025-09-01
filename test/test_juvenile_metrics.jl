"""
Tests for juvenile metrics
"""

using Test
using ReefMetrics

@testset "Relative Juveniles over Locations" begin
    n_tsteps, n_groups, n_sizes, n_locs = 2, 2, 4, 2
    relative_cover = zeros(Float64, n_tsteps, n_groups, n_sizes, n_locs)

    is_juvenile = [true, true, false, false]

    relative_cover[1, 1, :, 1] = [0.1, 0.1, 0.05, 0.05]
    relative_cover[1, 2, :, 1] = [0.05, 0.05, 0.2, 0.2]

    relative_cover[2, 1, :, 2] = [0.0, 0.05, 0.3, 0.3]
    relative_cover[2, 2, :, 2] = [0.2, 0.0, 0.05, 0.0]

    @testset "Normal Cases" begin
        rel_juv = relative_loc_juveniles(relative_cover, is_juvenile)
        @test size(rel_juv) == (n_tsteps, n_locs)

        @test rel_juv[1, 1] ≈ 0.3
        @test rel_juv[1, 2] == 0.0
        @test rel_juv[2, 1] == 0.0
        @test rel_juv[2, 2] == 0.25
    end

    @testset "Edge Cases" begin
        # Test with no juvenile classes
        no_juveniles = [false, false, false, false]
        @test all(relative_loc_juveniles(relative_cover, no_juveniles) .== 0.0)

        # Test with all classes as juveniles
        all_juveniles = [true, true, true, true]
        expected_sum = dropdims(sum(relative_cover, dims=(2, 3)), dims=(2, 3))
        @test relative_loc_juveniles(relative_cover, all_juveniles) ≈ expected_sum
    end
end

@testset "Absolute Juveniles" begin
    n_tsteps, n_groups, n_sizes, n_locs = 2, 2, 4, 2
    relative_cover = zeros(Float64, n_tsteps, n_groups, n_sizes, n_locs)
    location_area = [100.0, 200.0]

    is_juvenile = [true, true, false, false]

    relative_cover[1, 1, :, 1] = [0.1, 0.1, 0.05, 0.05]
    relative_cover[1, 2, :, 1] = [0.05, 0.05, 0.2, 0.2]

    relative_cover[2, 1, :, 2] = [0.0, 0.05, 0.3, 0.3]
    relative_cover[2, 2, :, 2] = [0.2, 0.0, 0.05, 0.0]

    @testset "Noraml Cases" begin
        abs_juv = absolute_loc_juveniles(relative_cover, is_juvenile, location_area)
        @test size(abs_juv) == (n_tsteps, n_locs)

        @test abs_juv[1, 1] ≈ 30.0
        @test abs_juv[1, 2] ≈ 0.0
        @test abs_juv[2, 1] ≈ 0.0
        @test abs_juv[2, 2] ≈ 50.0
    end

    @testset "Edge Cases" begin
        no_juveniles = [false, false, false, false]
        @test all(absolute_loc_juveniles(
            relative_cover, no_juveniles, location_area
        ) .== 0.0)

        all_juveniles = [true, true, true, true]
        abs_juv = absolute_loc_juveniles(relative_cover, all_juveniles, location_area)

        @test abs_juv[1, 1] ≈ 80.0
        @test abs_juv[1, 2] ≈ 0.0
        @test abs_juv[2, 1] ≈ 0.0
        @test abs_juv[2, 2] ≈ 180.0
    end
end

@testset "Juvenile Indicator" begin
    n_tsteps, n_groups, n_sizes, n_locs = 2, 2, 4, 2
    relative_cover = zeros(Float64, n_tsteps, n_groups, n_sizes, n_locs)
    location_area = [100.0, 200.0]

    is_juvenile = [true, true, false, false]

    relative_cover[1, 1, :, 1] = [0.1, 0.1, 0.05, 0.05]
    relative_cover[1, 2, :, 1] = [0.05, 0.05, 0.2, 0.2]

    relative_cover[2, 1, :, 2] = [0.0, 0.05, 0.3, 0.3]
    relative_cover[2, 2, :, 2] = [0.2, 0.0, 0.05, 0.0]

    @testset "Normal Cases" begin
        max_juv_colony_area = 0.1
        max_juv_density = 3.0
        juv_ind = juvenile_indicator(
            relative_cover, is_juvenile, location_area, max_juv_colony_area, max_juv_density
        )

        @test juv_ind[1, 1] ≈ 30.0 / 30.0
        @test juv_ind[1, 2] ≈ 0.0
        @test juv_ind[2, 1] ≈ 0.0
        @test juv_ind[2, 2] ≈ 50.0 / 60.0

        max_juv_colony_area = 0.2
        max_juv_density = 3.0
        juv_ind = juvenile_indicator(
            relative_cover, is_juvenile, location_area, max_juv_colony_area, max_juv_density
        )

        @test juv_ind[1, 1] ≈ 30.0 / 60.0
        @test juv_ind[1, 2] ≈ 0.0
        @test juv_ind[2, 1] ≈ 0.0
        @test juv_ind[2, 2] ≈ 50.0 / 120.0
    end

    @testset "Edge Cases" begin
        no_juveniles = [false, false, false, false]
        max_juv_colony_area = 0.2
        max_juv_density = 3.0
        @test all(
            juvenile_indicator(
                relative_cover, no_juveniles, location_area, max_juv_colony_area,
                max_juv_density
            ) .== 0.0
        )
    end
end

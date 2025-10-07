"""
Tests for cover metrics
"""

using Test
using ADRIAIndicators: relative_cover, relative_taxa_cover, relative_loc_taxa_cover, ltmp_cover, ltmp_taxa_cover, ltmp_loc_taxa_cover

@testset "Relative Cover" begin
    n_tsteps, n_groups, n_sizes, n_locs = 2, 2, 3, 2
    relative_cover_data = zeros(Float64, n_tsteps, n_groups, n_sizes, n_locs)

    # Timestep 1, Location 1
    relative_cover_data[1, 1, :, 1] = [0.1, 0.1, 0.05]  # Group 1 -> sum = 0.25
    relative_cover_data[1, 2, :, 1] = [0.05, 0.05, 0.2] # Group 2 -> sum = 0.3
    # Total for loc 1, time 1 = 0.25 + 0.3 = 0.55

    # Timestep 2, Location 2
    relative_cover_data[2, 1, :, 2] = [0.0, 0.05, 0.3]   # Group 1 -> sum = 0.35
    relative_cover_data[2, 2, :, 2] = [0.2, 0.0, 0.05]   # Group 2 -> sum = 0.25
    # Total for loc 2, time 2 = 0.35 + 0.25 = 0.6

    @testset "Normal Cases" begin
        rc = relative_cover(relative_cover_data)
        @test size(rc) == (n_tsteps, n_locs)

        @test rc[1, 1] ≈ 0.55
        @test rc[1, 2] == 0.0
        @test rc[2, 1] == 0.0
        @test rc[2, 2] ≈ 0.6
    end

    @testset "Edge Cases" begin
        # All zeros
        @test all(relative_cover(zeros(Float64, 2, 2, 2, 2)) .== 0.0)
    end
end

@testset "Relative Taxa Cover" begin
    n_tsteps, n_groups, n_sizes, n_locs = 2, 2, 3, 2
    relative_cover_data = zeros(Float64, n_tsteps, n_groups, n_sizes, n_locs)
    location_area = [100.0, 200.0]
    total_area = sum(location_area)

    # Timestep 1, Location 1
    relative_cover_data[1, 1, :, 1] = [0.1, 0.1, 0.05]  # Group 1 -> sum = 0.25
    relative_cover_data[1, 2, :, 1] = [0.05, 0.05, 0.2] # Group 2 -> sum = 0.3

    # Timestep 2, Location 2
    relative_cover_data[2, 1, :, 2] = [0.0, 0.05, 0.3]   # Group 1 -> sum = 0.35
    relative_cover_data[2, 2, :, 2] = [0.2, 0.0, 0.05]   # Group 2 -> sum = 0.25

    @testset "Normal Cases" begin
        rtc = relative_taxa_cover(relative_cover_data, location_area)
        @test size(rtc) == (n_tsteps, n_groups)

        # Expected values are (sum of absolute taxa cover) / (total area)
        # T1, G1: (0.25 * 100) / 300
        @test rtc[1, 1] ≈ (0.25 * 100) / total_area
        # T1, G2: (0.3 * 100) / 300
        @test rtc[1, 2] ≈ (0.3 * 100) / total_area

        # T2, G1: (0.35 * 200) / 300
        @test rtc[2, 1] ≈ (0.35 * 200) / total_area
        # T2, G2: (0.25 * 200) / 300
        @test rtc[2, 2] ≈ (0.25 * 200) / total_area
    end

    @testset "DimensionMismatch" begin
        @test_throws DimensionMismatch relative_taxa_cover(relative_cover_data, [100.0])
    end
end

@testset "Relative Location Taxa Cover" begin
    n_tsteps, n_groups, n_sizes, n_locs = 2, 2, 3, 2
    relative_cover_data = zeros(Float64, n_tsteps, n_groups, n_sizes, n_locs)

    # Timestep 1, Location 1
    relative_cover_data[1, 1, :, 1] = [0.1, 0.1, 0.05]  # Group 1 -> sum = 0.25
    relative_cover_data[1, 2, :, 1] = [0.05, 0.05, 0.2] # Group 2 -> sum = 0.3

    # Timestep 2, Location 2
    relative_cover_data[2, 1, :, 2] = [0.0, 0.05, 0.3]   # Group 1 -> sum = 0.35
    relative_cover_data[2, 2, :, 2] = [0.2, 0.0, 0.05]   # Group 2 -> sum = 0.25

    @testset "Normal Cases" begin
        rltc = relative_loc_taxa_cover(relative_cover_data)
        @test size(rltc) == (n_tsteps, n_groups, n_locs)

        @test rltc[1, 1, 1] ≈ 0.25
        @test rltc[1, 2, 1] ≈ 0.3
        @test all(rltc[1, :, 2] .== 0.0)

        @test all(rltc[2, :, 1] .== 0.0)
        @test rltc[2, 1, 2] ≈ 0.35
        @test rltc[2, 2, 2] ≈ 0.25
    end
end

@testset "LTMP Cover Metrics" begin
    n_tsteps, n_groups, n_sizes, n_locs = 2, 2, 3, 2
    relative_cover_data = zeros(Float64, n_tsteps, n_groups, n_sizes, n_locs)

    # Timestep 1, Location 1
    relative_cover_data[1, 1, :, 1] = [0.1, 0.1, 0.05]  # Group 1 -> sum = 0.25
    relative_cover_data[1, 2, :, 1] = [0.05, 0.05, 0.2] # Group 2 -> sum = 0.3
    # Total relative cover for loc 1, time 1 = 0.55

    # Timestep 2, Location 2
    relative_cover_data[2, 1, :, 2] = [0.0, 0.05, 0.3]   # Group 1 -> sum = 0.35
    relative_cover_data[2, 2, :, 2] = [0.2, 0.0, 0.05]   # Group 2 -> sum = 0.25
    # Total relative cover for loc 2, time 2 = 0.6

    habitable_area = [100.0, 200.0]
    reef_area = [150.0, 400.0]

    @testset "LTMP Cover" begin
        ltmp = ltmp_cover(relative_cover_data, habitable_area, reef_area)
        @test size(ltmp) == (n_tsteps, n_locs)

        # Expected: (relative_cover * habitable_area) / reef_area
        @test ltmp[1, 1] ≈ (0.55 * 100.0) / 150.0
        @test ltmp[1, 2] == 0.0
        @test ltmp[2, 1] == 0.0
        @test ltmp[2, 2] ≈ (0.6 * 200.0) / 400.0
    end

    @testset "LTMP Taxa Cover" begin
        total_reef_area = sum(reef_area)
        ltmp_tc = ltmp_taxa_cover(relative_cover_data, habitable_area, reef_area)
        @test size(ltmp_tc) == (n_tsteps, n_groups)

        # Expected: (sum of absolute taxa cover) / (total reef area)
        # T1, G1: (0.25 * 100) / 550
        @test ltmp_tc[1, 1] ≈ (0.25 * 100.0) / total_reef_area
        # T1, G2: (0.3 * 100) / 550
        @test ltmp_tc[1, 2] ≈ (0.3 * 100.0) / total_reef_area

        # T2, G1: (0.35 * 200) / 550
        @test ltmp_tc[2, 1] ≈ (0.35 * 200.0) / total_reef_area
        # T2, G2: (0.25 * 200) / 550
        @test ltmp_tc[2, 2] ≈ (0.25 * 200.0) / total_reef_area
    end

    @testset "LTMP Location Taxa Cover" begin
        ltmp_ltc = ltmp_loc_taxa_cover(relative_cover_data, habitable_area, reef_area)
        @test size(ltmp_ltc) == (n_tsteps, n_groups, n_locs)

        # Expected: relative_loc_taxa_cover * (habitable_area / reef_area)
        @test ltmp_ltc[1, 1, 1] ≈ 0.25 * (100.0 / 150.0)
        @test ltmp_ltc[1, 2, 1] ≈ 0.3 * (100.0 / 150.0)
        @test all(ltmp_ltc[1, :, 2] .== 0.0)

        @test all(ltmp_ltc[2, :, 1] .== 0.0)
        @test ltmp_ltc[2, 1, 2] ≈ 0.35 * (200.0 / 400.0)
        @test ltmp_ltc[2, 2, 2] ≈ 0.25 * (200.0 / 400.0)
    end
end
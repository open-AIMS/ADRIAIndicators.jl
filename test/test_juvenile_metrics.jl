"""
Tests for juvenile metrics
"""

using Test
using ReefMetrics

@testset "Relative Juveniles" begin
    n_tsteps, n_groups, n_sizes, n_locs = 2, 2, 4, 2
    relative_cover = zeros(Float64, n_tsteps, n_groups, n_sizes, n_locs)
    location_area = [100.0, 200.0]
    total_area = sum(location_area)

    is_juvenile = [true, true, false, false]

    # Timestep 1: Cover only exists at location 1
    relative_cover[1, 1, :, 1] = [0.1, 0.1, 0.05, 0.05]
    relative_cover[1, 2, :, 1] = [0.05, 0.05, 0.2, 0.2]

    # Timestep 2: Cover only exists at location 2
    relative_cover[2, 1, :, 2] = [0.0, 0.05, 0.3, 0.3]
    relative_cover[2, 2, :, 2] = [0.2, 0.0, 0.05, 0.0]

    @testset "Normal Cases" begin
        rel_juv = relative_juveniles(relative_cover, is_juvenile, location_area)

        # Expected output is a vector of length n_tsteps
        @test size(rel_juv) == (n_tsteps,)

        # Calculation for Timestep 1:
        # Absolute juvenile cover at loc 1: ((0.1 + 0.1) + (0.05 + 0.05)) * 100.0 = 0.3 * 100.0 = 30.0
        # Absolute juvenile cover at loc 2: 0.0
        # Total absolute juvenile cover: 30.0
        # Relative to total area: 30.0 / 300.0 = 0.1
        @test rel_juv[1] ≈ 30.0 / total_area

        # Calculation for Timestep 2:
        # Absolute juvenile cover at loc 1: 0.0
        # Absolute juvenile cover at loc 2: ((0.0 + 0.05) + (0.2 + 0.0)) * 200.0 = 0.25 * 200.0 = 50.0
        # Total absolute juvenile cover: 50.0
        # Relative to total area: 50.0 / 300.0
        @test rel_juv[2] ≈ 50.0 / total_area
    end

    @testset "Edge Cases" begin
        # Test with no juvenile classes
        no_juveniles = [false, false, false, false]
        @test all(relative_juveniles(relative_cover, no_juveniles, location_area) .== 0.0)

        # Test with all classes as juveniles
        all_juveniles = [true, true, true, true]
        rel_juv_all = relative_juveniles(relative_cover, all_juveniles, location_area)

        # Calculation for Timestep 1 (all juveniles):
        # Absolute cover at loc 1: sum(relative_cover[1, :, :, 1]) * 100.0 = 0.8 * 100.0 = 80.0
        # Relative to total area: 80.0 / 300.0
        @test rel_juv_all[1] ≈ (sum(relative_cover[1, :, :, 1]) * location_area[1]) / total_area

        # Calculation for Timestep 2 (all juveniles):
        # Absolute cover at loc 2: sum(relative_cover[2, :, :, 2]) * 200.0 = 0.9 * 200.0 = 180.0
        # Relative to total area: 180.0 / 300.0
        @test rel_juv_all[2] ≈ (sum(relative_cover[2, :, :, 2]) * location_area[2]) / total_area

        # Test with zero area for all locations
        zero_area = [0.0, 0.0]
        # Division by zero should result in NaN
        @test all(isnan.(relative_juveniles(relative_cover, is_juvenile, zero_area)))
    end
end

@testset "Relative Taxa Juveniles" begin
    n_tsteps, n_groups, n_sizes, n_locs = 2, 2, 4, 2
    relative_cover = zeros(Float64, n_tsteps, n_groups, n_sizes, n_locs)
    location_area = [100.0, 200.0]
    total_area = sum(location_area)

    is_juvenile = [true, true, false, false]

    # Timestep 1: Cover is only in location 1
    # Group 1, Loc 1: Juvenile cover = 0.1 + 0.1 = 0.2. Absolute = 0.2 * 100 = 20
    relative_cover[1, 1, :, 1] = [0.1, 0.1, 0.05, 0.05]
    # Group 2, Loc 1: Juvenile cover = 0.05 + 0.05 = 0.1. Absolute = 0.1 * 100 = 10
    relative_cover[1, 2, :, 1] = [0.05, 0.05, 0.2, 0.2]

    # Timestep 2: Cover is only in location 2
    # Group 1, Loc 2: Juvenile cover = 0.0 + 0.05 = 0.05. Absolute = 0.05 * 200 = 10
    relative_cover[2, 1, :, 2] = [0.0, 0.05, 0.3, 0.3]
    # Group 2, Loc 2: Juvenile cover = 0.2 + 0.0 = 0.2. Absolute = 0.2 * 200 = 40
    relative_cover[2, 2, :, 2] = [0.2, 0.0, 0.05, 0.0]

    @testset "Normal Cases" begin
        rel_taxa_juv = relative_taxa_juveniles(relative_cover, is_juvenile, location_area)

        @test size(rel_taxa_juv) == (n_tsteps, n_groups)

        # Expected values are (sum of absolute juvenile cover) / (total area)
        # Timestep 1, Group 1: (0.2 * 100) / 300 = 20 / 300
        @test rel_taxa_juv[1, 1] ≈ 20.0 / total_area
        # Timestep 1, Group 2: (0.1 * 100) / 300 = 10 / 300
        @test rel_taxa_juv[1, 2] ≈ 10.0 / total_area

        # Timestep 2, Group 1: (0.05 * 200) / 300 = 10 / 300
        @test rel_taxa_juv[2, 1] ≈ 10.0 / total_area
        # Timestep 2, Group 2: (0.2 * 200) / 300 = 40 / 300
        @test rel_taxa_juv[2, 2] ≈ 40.0 / total_area
    end

    @testset "Edge Cases" begin
        # Test with no juvenile classes
        no_juveniles = [false, false, false, false]
        @test all(
            relative_taxa_juveniles(relative_cover, no_juveniles, location_area) .== 0.0
        )

        # Test with all classes as juveniles
        all_juveniles = [true, true, true, true]
        rel_taxa_all_juv = relative_taxa_juveniles(
            relative_cover, all_juveniles, location_area
        )

        # Calculate expected total relative cover per taxon (weighted by area)
        # T1, G1: (sum([0.1, 0.1, 0.05, 0.05])) * 100 / 300
        @test rel_taxa_all_juv[1, 1] ≈ 0.1
        # T1, G2: (sum([0.05, 0.05, 0.2, 0.2])) * 100 / 300
        @test rel_taxa_all_juv[1, 2] ≈ 5 / 30
        # T2, G1: (sum([0.0, 0.05, 0.3, 0.3])) * 200 / 300
        @test rel_taxa_all_juv[2, 1] ≈ 1.3 / 3
        # T2, G2: (sum([0.2, 0.0, 0.05, 0.0])) * 200 / 300
        @test rel_taxa_all_juv[2, 2] ≈ 5 / 30
    end
end

@testset "Relative Location Juveniles" begin
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

@testset "Relative Location Taxa Juveniles" begin
    n_tsteps, n_groups, n_sizes, n_locs = 2, 2, 4, 2
    relative_cover = zeros(Float64, n_tsteps, n_groups, n_sizes, n_locs)

    is_juvenile = [true, true, false, false]

    # Timestep 1, Location 1
    relative_cover[1, 1, :, 1] .= [0.1, 0.1, 0.05, 0.05]  # Group 1
    relative_cover[1, 2, :, 1] .= [0.05, 0.05, 0.2, 0.2]   # Group 2

    # Timestep 2, Location 2
    relative_cover[2, 1, :, 2] .= [0.0, 0.05, 0.3, 0.3]   # Group 1
    relative_cover[2, 2, :, 2] .= [0.2, 0.0, 0.05, 0.0]   # Group 2

    @testset "Normal Cases" begin
        rel_juv_taxa = relative_loc_taxa_juveniles(relative_cover, is_juvenile)

        # Expected shape is [timesteps, groups, locations]
        @test size(rel_juv_taxa) == (n_tsteps, n_groups, n_locs)

        @test rel_juv_taxa[1, 1, 1] ≈ 0.2
        @test rel_juv_taxa[1, 2, 1] ≈ 0.1
        @test all(rel_juv_taxa[1, :, 2] .== 0.0)
        @test all(rel_juv_taxa[2, :, 1] .== 0.0)
        @test rel_juv_taxa[2, 1, 2] ≈ 0.05
        @test rel_juv_taxa[2, 2, 2] ≈ 0.2
    end

    @testset "Edge Cases" begin
        # Test with no juvenile classes
        no_juveniles = [false, false, false, false]
        @test all(
            relative_loc_taxa_juveniles(relative_cover, no_juveniles) .== 0.0
        )

        # Test with all classes as juveniles
        all_juveniles = [true, true, true, true]
        rel_juv_all = relative_loc_taxa_juveniles(
            relative_cover, all_juveniles
        )

        @test rel_juv_all[1, 1, 1] ≈ 0.3
    end
end

@testset "Absolute Juveniles" begin
    n_tsteps, n_groups, n_sizes, n_locs = 2, 2, 4, 2
    relative_cover = zeros(Float64, n_tsteps, n_groups, n_sizes, n_locs)
    location_area = [100.0, 200.0]
    total_area = sum(location_area)

    is_juvenile = [true, true, false, false]

    # Timestep 1: Cover only exists at location 1
    relative_cover[1, 1, :, 1] = [0.1, 0.1, 0.05, 0.05]
    relative_cover[1, 2, :, 1] = [0.05, 0.05, 0.2, 0.2]

    # Timestep 2: Cover only exists at location 2
    relative_cover[2, 1, :, 2] = [0.0, 0.05, 0.3, 0.3]
    relative_cover[2, 2, :, 2] = [0.2, 0.0, 0.05, 0.0]

    @testset "Normal Cases" begin
        rel_juv = absolute_juveniles(relative_cover, is_juvenile, location_area)

        # Expected output is a vector of length n_tsteps
        @test size(rel_juv) == (n_tsteps,)

        @test rel_juv[1] ≈ 30.0
        @test rel_juv[2] ≈ 50.0
    end

    @testset "Edge Cases" begin
        # Test with no juvenile classes
        no_juveniles = [false, false, false, false]
        @test all(absolute_juveniles(relative_cover, no_juveniles, location_area) .== 0.0)

        # Test with all classes as juveniles
        all_juveniles = [true, true, true, true]
        rel_juv_all = absolute_juveniles(relative_cover, all_juveniles, location_area)

        # Calculation for Timestep 1 (all juveniles):
        # Absolute cover at loc 1: sum(relative_cover[1, :, :, 1]) * 100.0 = 0.8 * 100.0 = 80.0
        @test rel_juv_all[1] ≈ (sum(relative_cover[1, :, :, 1]) * location_area[1])

        # Calculation for Timestep 2 (all juveniles):
        # Absolute cover at loc 2: sum(relative_cover[2, :, :, 2]) * 200.0 = 0.9 * 200.0 = 180.0
        @test rel_juv_all[2] ≈ (sum(relative_cover[2, :, :, 2]) * location_area[2])

        # Test with zero area for all locations
        zero_area = [0.0, 0.0]
        @test all(absolute_juveniles(relative_cover, is_juvenile, zero_area) .== 0.0)
    end
end

@testset "Absolute Taxa Juveniles" begin
    n_tsteps, n_groups, n_sizes, n_locs = 2, 2, 4, 2
    relative_cover = zeros(Float64, n_tsteps, n_groups, n_sizes, n_locs)
    location_area = [100.0, 200.0]
    total_area = sum(location_area)

    is_juvenile = [true, true, false, false]

    # Timestep 1: Cover is only in location 1
    # Group 1, Loc 1: Juvenile cover = 0.1 + 0.1 = 0.2. Absolute = 0.2 * 100 = 20
    relative_cover[1, 1, :, 1] = [0.1, 0.1, 0.05, 0.05]
    # Group 2, Loc 1: Juvenile cover = 0.05 + 0.05 = 0.1. Absolute = 0.1 * 100 = 10
    relative_cover[1, 2, :, 1] = [0.05, 0.05, 0.2, 0.2]

    # Timestep 2: Cover is only in location 2
    # Group 1, Loc 2: Juvenile cover = 0.0 + 0.05 = 0.05. Absolute = 0.05 * 200 = 10
    relative_cover[2, 1, :, 2] = [0.0, 0.05, 0.3, 0.3]
    # Group 2, Loc 2: Juvenile cover = 0.2 + 0.0 = 0.2. Absolute = 0.2 * 200 = 40
    relative_cover[2, 2, :, 2] = [0.2, 0.0, 0.05, 0.0]

    @testset "Normal Cases" begin
        rel_taxa_juv = absolute_taxa_juveniles(relative_cover, is_juvenile, location_area)

        @test size(rel_taxa_juv) == (n_tsteps, n_groups)
        @test rel_taxa_juv[1, 1] ≈ 20.0
        @test rel_taxa_juv[1, 2] ≈ 10.0
        @test rel_taxa_juv[2, 1] ≈ 10.0
        @test rel_taxa_juv[2, 2] ≈ 40.0
    end

    @testset "Edge Cases" begin
        # Test with no juvenile classes
        no_juveniles = [false, false, false, false]
        @test all(
            absolute_taxa_juveniles(relative_cover, no_juveniles, location_area) .== 0.0
        )

        # Test with all classes as juveniles
        all_juveniles = [true, true, true, true]
        rel_taxa_all_juv = absolute_taxa_juveniles(
            relative_cover, all_juveniles, location_area
        )

        # Calculate expected total relative cover per taxon (weighted by area)
        # T1, G1: (sum([0.1, 0.1, 0.05, 0.05])) * 100 / 300
        @test rel_taxa_all_juv[1, 1] ≈ 30.0
        # T1, G2: (sum([0.05, 0.05, 0.2, 0.2])) * 100 / 300
        @test rel_taxa_all_juv[1, 2] ≈ 50.0
        # T2, G1: (sum([0.0, 0.05, 0.3, 0.3])) * 200 / 300
        @test rel_taxa_all_juv[2, 1] ≈ 130.0
        # T2, G2: (sum([0.2, 0.0, 0.05, 0.0])) * 200 / 300
        @test rel_taxa_all_juv[2, 2] ≈ 50.0
    end
end

@testset "Absolute Location Juveniles" begin
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

@testset "Absolute Location Taxa Juveniles" begin
    n_tsteps, n_groups, n_sizes, n_locs = 2, 2, 4, 2
    relative_cover = zeros(Float64, n_tsteps, n_groups, n_sizes, n_locs)
    location_area = [100.0, 200.0]

    is_juvenile = [true, true, false, false]

    # Timestep 1, Location 1
    relative_cover[1, 1, :, 1] .= [0.1, 0.1, 0.05, 0.05]  # Group 1
    relative_cover[1, 2, :, 1] .= [0.05, 0.05, 0.2, 0.2]   # Group 2

    # Timestep 2, Location 2
    relative_cover[2, 1, :, 2] .= [0.0, 0.05, 0.3, 0.3]   # Group 1
    relative_cover[2, 2, :, 2] .= [0.2, 0.0, 0.05, 0.0]   # Group 2

    @testset "Normal Cases" begin
        rel_juv_taxa = absolute_loc_taxa_juveniles(relative_cover, is_juvenile, location_area)

        # Expected shape is [timesteps, groups, locations]
        @test size(rel_juv_taxa) == (n_tsteps, n_groups, n_locs)

        @test rel_juv_taxa[1, 1, 1] ≈ 0.2 * 100.0
        @test rel_juv_taxa[1, 2, 1] ≈ 0.1 * 100.0
        @test all(rel_juv_taxa[1, :, 2] .== 0.0)
        @test all(rel_juv_taxa[2, :, 1] .== 0.0)
        @test rel_juv_taxa[2, 1, 2] ≈ 0.05 * 200.0
        @test rel_juv_taxa[2, 2, 2] ≈ 0.2 * 200.0
    end

    @testset "Edge Cases" begin
        # Test with no juvenile classes
        no_juveniles = [false, false, false, false]
        @test all(
            absolute_loc_taxa_juveniles(relative_cover, no_juveniles, location_area) .== 0.0
        )

        # Test with all classes as juveniles
        all_juveniles = [true, true, true, true]
        rel_juv_all = absolute_loc_taxa_juveniles(
            relative_cover, all_juveniles, location_area
        )

        @test rel_juv_all[1, 1, 1] ≈ 0.3 * 100.0
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

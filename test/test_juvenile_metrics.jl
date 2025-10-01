"""
Tests for juvenile metrics
"""

using Test
using ADRIAIndicators: relative_juveniles, relative_loc_taxa_juveniles,
    relative_taxa_juveniles

using ADRIAIndicators: absolute_juveniles, absolute_loc_taxa_juveniles,
    absolute_taxa_juveniles

using ADRIAIndicators: juvenile_indicator

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

@testset "Relative Juveniles" begin
    n_tsteps, n_groups, n_sizes, n_locs = 2, 2, 4, 2
    relative_cover = zeros(Float64, n_tsteps, n_groups, n_sizes, n_locs)

    is_juvenile = [true, true, false, false]

    relative_cover[1, 1, :, 1] = [0.1, 0.1, 0.05, 0.05]
    relative_cover[1, 2, :, 1] = [0.05, 0.05, 0.2, 0.2]

    relative_cover[2, 1, :, 2] = [0.0, 0.05, 0.3, 0.3]
    relative_cover[2, 2, :, 2] = [0.2, 0.0, 0.05, 0.0]

    @testset "Normal Cases" begin
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
        expected_sum = dropdims(sum(relative_cover, dims=(2, 3)), dims=(2, 3))
        @test relative_juveniles(relative_cover, all_juveniles) ≈ expected_sum
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
        abs_juv = absolute_juveniles(relative_cover, is_juvenile, location_area)
        @test size(abs_juv) == (n_tsteps, n_locs)

        @test abs_juv[1, 1] ≈ 30.0
        @test abs_juv[1, 2] ≈ 0.0
        @test abs_juv[2, 1] ≈ 0.0
        @test abs_juv[2, 2] ≈ 50.0
    end

    @testset "Edge Cases" begin
        no_juveniles = [false, false, false, false]
        @test all(absolute_juveniles(
            relative_cover, no_juveniles, location_area
        ) .== 0.0)

        all_juveniles = [true, true, true, true]
        abs_juv = absolute_juveniles(relative_cover, all_juveniles, location_area)

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
    habitable_area = [100.0, 200.0]

    is_juvenile = [true, true, false, false]

    relative_cover[1, 1, :, 1] = [0.1, 0.1, 0.05, 0.05]
    relative_cover[1, 2, :, 1] = [0.05, 0.05, 0.2, 0.2]

    relative_cover[2, 1, :, 2] = [0.0, 0.05, 0.3, 0.3]
    relative_cover[2, 2, :, 2] = [0.2, 0.0, 0.05, 0.0]

    mean_colony_diameters = zeros(Float64, n_groups, n_sizes)
    mean_colony_diameters[1, :] = [0.1, sqrt(0.4 / π), 0.5, 0.6]
    mean_colony_diameters[2, :] = [0.1, 0.1, 0.5, 0.6]

    @testset "Normal Cases" begin
        max_juv_density = 3.0
        juv_ind = juvenile_indicator(
            relative_cover, is_juvenile, habitable_area, mean_colony_diameters, max_juv_density
        )

        # With mean_colony_diameters, max juvenile diameter is sqrt(0.4 / π)
        # This makes max_col_area = (π / 4) * (sqrt(0.4 / π))^2 = 0.1
        # This is to match the previous test's max_juv_colony_area = 0.1
        # Denominator for loc 1: 0.1 * 3.0 * 100.0 = 30.0
        # Denominator for loc 2: 0.1 * 3.0 * 200.0 = 60.0
        @test juv_ind[1, 1] ≈ 30.0 / 30.0
        @test juv_ind[1, 2] ≈ 0.0
        @test juv_ind[2, 1] ≈ 0.0
        @test juv_ind[2, 2] ≈ 50.0 / 60.0

        # Another test with different diameters
        mean_colony_diameters[1, 2] = sqrt(0.8 / π)
        # Now max juvenile diameter is sqrt(0.8 / π)
        # max_col_area is now 0.2
        juv_ind = juvenile_indicator(
            relative_cover, is_juvenile, habitable_area, mean_colony_diameters, max_juv_density
        )

        # Denominator for loc 1: 0.2 * 3.0 * 100.0 = 60.0
        # Denominator for loc 2: 0.2 * 3.0 * 200.0 = 120.0
        @test juv_ind[1, 1] ≈ 30.0 / 60.0
        @test juv_ind[1, 2] ≈ 0.0
        @test juv_ind[2, 1] ≈ 0.0
        @test juv_ind[2, 2] ≈ 50.0 / 120.0
    end

    @testset "Edge Cases" begin
        no_juveniles = [false, false, false, false]
        max_juv_density = 3.0
        @test all(
            juvenile_indicator(
                relative_cover, no_juveniles, habitable_area, mean_colony_diameters,
                max_juv_density
            ) .== 0.0
        )
    end
end

@testset "Relative Juveniles 5D" begin
    n_tsteps, n_groups, n_sizes, n_locs, n_scenarios = 2, 2, 4, 2, 2
    relative_cover = zeros(Float64, n_tsteps, n_groups, n_sizes, n_locs, n_scenarios)

    is_juvenile = [true, true, false, false]

    # Scenario 1
    relative_cover[1, 1, :, 1, 1] = [0.1, 0.1, 0.05, 0.05]
    relative_cover[1, 2, :, 1, 1] = [0.05, 0.05, 0.2, 0.2]
    relative_cover[2, 1, :, 2, 1] = [0.0, 0.05, 0.3, 0.3]
    relative_cover[2, 2, :, 2, 1] = [0.2, 0.0, 0.05, 0.0]

    # Scenario 2 - different values
    relative_cover[1, 1, :, 1, 2] = [0.2, 0.2, 0.1, 0.1]
    relative_cover[1, 2, :, 1, 2] = [0.1, 0.1, 0.3, 0.3]
    relative_cover[2, 1, :, 2, 2] = [0.1, 0.1, 0.4, 0.4]
    relative_cover[2, 2, :, 2, 2] = [0.3, 0.1, 0.1, 0.1]

    rel_juv_5d = relative_juveniles(relative_cover, is_juvenile)
    @test size(rel_juv_5d) == (n_tsteps, n_locs, n_scenarios)

    # Check scenario 1 (same as 4D test)
    @test rel_juv_5d[1, 1, 1] ≈ 0.3
    @test rel_juv_5d[1, 2, 1] == 0.0
    @test rel_juv_5d[2, 1, 1] == 0.0
    @test rel_juv_5d[2, 2, 1] == 0.25

    # Check scenario 2
    @test rel_juv_5d[1, 1, 2] ≈ 0.6
    @test rel_juv_5d[1, 2, 2] == 0.0
    @test rel_juv_5d[2, 1, 2] == 0.0
    @test rel_juv_5d[2, 2, 2] ≈ 0.6
end

@testset "Absolute Juveniles 5D" begin
    n_tsteps, n_groups, n_sizes, n_locs, n_scenarios = 2, 2, 4, 2, 2
    relative_cover = zeros(Float64, n_tsteps, n_groups, n_sizes, n_locs, n_scenarios)
    location_area = [100.0, 200.0]

    is_juvenile = [true, true, false, false]

    # Scenario 1
    relative_cover[1, 1, :, 1, 1] = [0.1, 0.1, 0.05, 0.05]
    relative_cover[1, 2, :, 1, 1] = [0.05, 0.05, 0.2, 0.2]
    relative_cover[2, 1, :, 2, 1] = [0.0, 0.05, 0.3, 0.3]
    relative_cover[2, 2, :, 2, 1] = [0.2, 0.0, 0.05, 0.0]

    # Scenario 2 - different values
    relative_cover[1, 1, :, 1, 2] = [0.2, 0.2, 0.1, 0.1]
    relative_cover[1, 2, :, 1, 2] = [0.1, 0.1, 0.3, 0.3]
    relative_cover[2, 1, :, 2, 2] = [0.1, 0.1, 0.4, 0.4]
    relative_cover[2, 2, :, 2, 2] = [0.3, 0.1, 0.1, 0.1]

    abs_juv_5d = absolute_juveniles(relative_cover, is_juvenile, location_area)
    @test size(abs_juv_5d) == (n_tsteps, n_locs, n_scenarios)

    # Check scenario 1 (same as 4D test)
    @test abs_juv_5d[1, 1, 1] ≈ 30.0
    @test abs_juv_5d[1, 2, 1] ≈ 0.0
    @test abs_juv_5d[2, 1, 1] ≈ 0.0
    @test abs_juv_5d[2, 2, 1] ≈ 50.0

    # Check scenario 2
    @test abs_juv_5d[1, 1, 2] ≈ 60.0
    @test abs_juv_5d[1, 2, 2] == 0.0
    @test abs_juv_5d[2, 1, 2] == 0.0
    @test abs_juv_5d[2, 2, 2] ≈ 120.0
end

@testset "Juvenile Indicator 5D" begin
    n_tsteps, n_groups, n_sizes, n_locs, n_scenarios = 2, 2, 4, 2, 2
    relative_cover = zeros(Float64, n_tsteps, n_groups, n_sizes, n_locs, n_scenarios)
    habitable_area = [100.0, 200.0]

    is_juvenile = [true, true, false, false]

    # Scenario 1
    relative_cover[1, 1, :, 1, 1] = [0.1, 0.1, 0.05, 0.05]
    relative_cover[1, 2, :, 1, 1] = [0.05, 0.05, 0.2, 0.2]
    relative_cover[2, 1, :, 2, 1] = [0.0, 0.05, 0.3, 0.3]
    relative_cover[2, 2, :, 2, 1] = [0.2, 0.0, 0.05, 0.0]

    # Scenario 2 - different values
    relative_cover[1, 1, :, 1, 2] = [0.2, 0.2, 0.1, 0.1]
    relative_cover[1, 2, :, 1, 2] = [0.1, 0.1, 0.3, 0.3]
    relative_cover[2, 1, :, 2, 2] = [0.1, 0.1, 0.4, 0.4]
    relative_cover[2, 2, :, 2, 2] = [0.3, 0.1, 0.1, 0.1]

    mean_colony_diameters = zeros(Float64, n_groups, n_sizes)
    mean_colony_diameters[1, :] = [0.1, sqrt(0.4 / π), 0.5, 0.6]
    mean_colony_diameters[2, :] = [0.1, 0.1, 0.5, 0.6]
    
    max_juv_density = 3.0

    juv_ind_5d = juvenile_indicator(
        relative_cover, is_juvenile, habitable_area, mean_colony_diameters, max_juv_density
    )

    @test size(juv_ind_5d) == (n_tsteps, n_locs, n_scenarios)

    # Check scenario 1 (same as 4D test)
    @test juv_ind_5d[1, 1, 1] ≈ 1.0
    @test juv_ind_5d[1, 2, 1] ≈ 0.0
    @test juv_ind_5d[2, 1, 1] ≈ 0.0
    @test juv_ind_5d[2, 2, 1] ≈ 50.0 / 60.0

    # Check scenario 2
    # abs_juv_5d[1, 1, 2] was 60.0
    # Denominator for loc 1: 0.1 * 3.0 * 100.0 = 30.0
    @test juv_ind_5d[1, 1, 2] ≈ 60.0 / 30.0
    @test juv_ind_5d[1, 2, 2] == 0.0
    @test juv_ind_5d[2, 1, 2] == 0.0
    # abs_juv_5d[2, 2, 2] was 120.0
    # Denominator for loc 2: 0.1 * 3.0 * 200.0 = 60.0
    @test juv_ind_5d[2, 2, 2] ≈ 120.0 / 60.0
end

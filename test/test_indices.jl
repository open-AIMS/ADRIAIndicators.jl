"""
Tests for reef indices.
"""

using Test
using ReefMetrics

@testset "Reef Condition Index" begin

    @testset "Very Good Condition" begin
        # All metrics are high enough to pass all thresholds
        rc = fill(0.5, 1, 1)
        sv = fill(0.5, 1, 1)
        juv = fill(0.4, 1, 1)
        cots = fill(0.0, 1, 1)    # High health contribution (cots_comp = 1.0)
        rubble = fill(0.1, 1, 1)  # High health contribution (rubble_comp = 0.9)

        expected_rci = [0.9]
        actual_rci = reef_condition_index(rc, sv, juv, cots, rubble)

        @test actual_rci ≈ expected_rci
    end

    @testset "Good Condition" begin
        # Metrics are set to fail the "Very Good" thresholds but pass "Good" and lower
        rc = fill(0.40, 1, 1)
        sv = fill(0.40, 1, 1)
        juv = fill(0.32, 1, 1)
        cots = fill(0.10, 1, 1)
        rubble = fill(0.18, 1, 1)

        expected_rci = [0.7]
        actual_rci = reef_condition_index(rc, sv, juv, cots, rubble)

        @test actual_rci ≈ expected_rci
    end

    @testset "Fair Condition" begin
        # Metrics are set to fail "Very Good" and "Good" but pass "Fair" and lower
        rc = fill(0.30, 1, 1)
        sv = fill(0.32, 1, 1)
        juv = fill(0.28, 1, 1)
        cots = fill(0.20, 1, 1)
        rubble = fill(0.22, 1, 1)

        expected_rci = [0.5]
        actual_rci = reef_condition_index(rc, sv, juv, cots, rubble)

        @test actual_rci ≈ expected_rci
    end

    @testset "Poor Condition" begin
        # Metrics are set to only pass the "Poor" threshold
        rc = fill(0.20, 1, 1)
        sv = fill(0.28, 1, 1)
        juv = fill(0.22, 1, 1)
        cots = fill(0.28, 1, 1)
        rubble = fill(0.28, 1, 1)

        expected_rci = [0.3]
        actual_rci = reef_condition_index(rc, sv, juv, cots, rubble)

        @test actual_rci ≈ expected_rci
    end

    @testset "Very Poor Condition" begin
        # Metrics are so low they fail all thresholds
        rc = fill(0.0, 1, 1)
        sv = fill(0.0, 1, 1)
        juv = fill(0.0, 1, 1)
        cots = fill(0.5, 1, 1)
        rubble = fill(0.9, 1, 1)

        expected_rci = [0.1]
        actual_rci = reef_condition_index(rc, sv, juv, cots, rubble)

        @test actual_rci ≈ expected_rci
    end

    @testset "Combined Matrix Test" begin
        # Each column represents a different condition category from VG to VP
        rc =     [0.5 0.40 0.30 0.20 0.0]
        sv =     [0.5 0.40 0.32 0.28 0.0]
        juv =    [0.4 0.32 0.28 0.22 0.0]
        cots =   [0.0 0.10 0.20 0.28 0.5]
        rubble = [0.1 0.18 0.22 0.28 0.9]

        expected_rci = [0.9 0.7 0.5 0.3 0.1]
        actual_rci = reef_condition_index(rc, sv, juv, cots, rubble)

        @test actual_rci ≈ expected_rci
    end
end

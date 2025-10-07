"""
Tests for reef indices.
"""

using Test
using ADRIAIndicators: reef_condition_index, reef_tourism_index,
    reef_biodiversity_condition_index, reef_fish_index

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

@testset "Reef Tourism Index" begin

    @testset "Core Calculation and Rounding" begin
        # Inputs chosen to produce a value that requires rounding.
        # Raw value = 0.47947 + 0.1 * (0.7678 + 0.2945 + 0.8371 + 0.2822 + 0.7764)
        #           = 0.47947 + 0.1 * 2.958 = 0.47947 + 0.2958 = 0.77527
        # Clamped value = 0.77527
        # Rounded value = 0.78
        rc =     fill(0.1, 1, 1)
        sv =     fill(0.1, 1, 1)
        juv =    fill(0.1, 1, 1)
        cots =   fill(0.1, 1, 1)
        rubble = fill(0.1, 1, 1)

        expected_rti = [0.78]

        # Test the wrapper function
        actual_rti_wrapper = reef_tourism_index(rc, sv, juv, cots, rubble)
        @test actual_rti_wrapper ≈ expected_rti
    end

    @testset "Dimension Mismatch" begin
        rc =     zeros(1, 2)
        sv =     zeros(1, 2)
        juv =    zeros(1, 2)
        cots =   zeros(1, 2)
        rubble = zeros(1, 1) # Mismatched dimension

        @test_throws DimensionMismatch reef_tourism_index(rc, sv, juv, cots, rubble)
    end

    @testset "Combined Matrix Test" begin
        # Each column represents a different test case:
        # 1. Standard/Rounding Case
        # 2. Upper Clamp Case
        # 3. Lower Clamp Case (with negative inputs)
        # 4. Zero Input Case
        rc =     [0.1 1.0 -1.0 0.0]
        sv =     [0.1 1.0 -1.0 0.0]
        juv =    [0.1 1.0 -1.0 0.0]
        cots =   [0.1 1.0 -1.0 0.0]
        rubble = [0.1 1.0 -1.0 0.0]

        expected_rti = [0.78 0.9 0.1 0.48]
        actual_rti = reef_tourism_index(rc, sv, juv, cots, rubble)

        @test actual_rti ≈ expected_rti
    end
end

@testset "Reef Condition Index (4 inputs)" begin

    @testset "Very Good Condition" begin
        # All metrics are high enough to pass all thresholds
        rc = fill(0.5, 1, 1)
        ce = fill(0.5, 1, 1)
        sv = fill(0.5, 1, 1)
        juv = fill(0.4, 1, 1)

        expected_rci = [0.9]
        actual_rci = reef_condition_index(rc, ce, sv, juv)

        @test actual_rci ≈ expected_rci
    end

    @testset "Good Condition" begin
        # Metrics are set to fail the "Very Good" thresholds but pass "Good" and lower
        rc = fill(0.40, 1, 1)
        ce = fill(0.40, 1, 1)
        sv = fill(0.32, 1, 1)
        juv = fill(0.28, 1, 1)

        expected_rci = [0.7]
        actual_rci = reef_condition_index(rc, ce, sv, juv)

        @test actual_rci ≈ expected_rci
    end

    @testset "Fair Condition" begin
        # Metrics are set to fail "Very Good" and "Good" but pass "Fair" and lower
        rc = fill(0.30, 1, 1)
        ce = fill(0.28, 1, 1)
        sv = fill(0.32, 1, 1)
        juv = fill(0.28, 1, 1)

        expected_rci = [0.5]
        actual_rci = reef_condition_index(rc, ce, sv, juv)

        @test actual_rci ≈ expected_rci
    end

    @testset "Poor Condition" begin
        # Metrics are set to only pass the "Poor" threshold
        rc = fill(0.20, 1, 1)
        ce = fill(0.26, 1, 1)
        sv = fill(0.31, 1, 1)
        juv = fill(0.20, 1, 1)

        expected_rci = [0.3]
        actual_rci = reef_condition_index(rc, ce, sv, juv)

        @test actual_rci ≈ expected_rci
    end

    @testset "Very Poor Condition" begin
        # Metrics are so low they fail all thresholds
        rc = fill(0.0, 1, 1)
        ce = fill(0.0, 1, 1)
        sv = fill(0.0, 1, 1)
        juv = fill(0.0, 1, 1)

        expected_rci = [0.1]
        actual_rci = reef_condition_index(rc, ce, sv, juv)

        @test actual_rci ≈ expected_rci
    end

    @testset "Combined Matrix Test" begin
        # Each column represents a different condition category from VG to VP
        rc =  [0.5 0.40 0.30 0.20 0.0]
        ce =  [0.5 0.40 0.28 0.26 0.0]
        sv =  [0.5 0.32 0.32 0.31 0.0]
        juv = [0.4 0.28 0.28 0.20 0.0]

        expected_rci = [0.9 0.7 0.5 0.3 0.1]
        actual_rci = reef_condition_index(rc, ce, sv, juv)

        @test actual_rci ≈ expected_rci
    end
end


@testset "Reef Tourism Index (4 inputs)" begin

    @testset "Core Calculation and Rounding" begin
        # Inputs chosen to produce a value that requires rounding.
        # Raw value = 0.47947 + 0.1 * (0.12764 + 0.31946 + 0.11676 - 0.0036065)
        #           = 0.47947 + 0.1 * 0.5602535 = 0.47947 + 0.05602535 = 0.53549535
        # Clamped value = 0.53549535
        # Rounded value = 0.54
        rc =     fill(0.1, 1, 1)
        ce =     fill(0.1, 1, 1)
        sv =     fill(0.1, 1, 1)
        juv =    fill(0.1, 1, 1)

        expected_rti = [0.54]

        actual_rti_wrapper = reef_tourism_index(rc, ce, sv, juv)
        @test actual_rti_wrapper ≈ expected_rti
    end

    @testset "Dimension Mismatch" begin
        rc =     zeros(1, 2)
        ce =     zeros(1, 2)
        sv =     zeros(1, 2)
        juv =    zeros(1, 1) # Mismatched dimension

        @test_throws DimensionMismatch reef_tourism_index(rc, ce, sv, juv)
    end

    @testset "Combined Matrix Test" begin
        # Each column represents a different test case:
        # 1. Standard/Rounding Case
        # 2. Upper Clamp Case
        # 3. Lower Clamp Case (with negative inputs)
        # 4. Zero Input Case
        rc =     [0.1 1.0 -1.0 0.0]
        ce =     [0.1 1.0 -1.0 0.0]
        sv =     [0.1 1.0 -1.0 0.0]
        juv =    [0.1 1.0 -1.0 0.0]

        # Expected values calculated from the linear model coefficients
        # Case 1: 0.54 (as above)
        # Case 2 (Upper Clamp): Raw = 0.47947 + (0.12764 + 0.31946 + 0.11676 - 0.0036065) = 1.0397235 -> 0.9
        # Case 3 (Lower Clamp): Raw = 0.47947 - 0.5602535 = -0.0807835 -> 0.1
        # Case 4 (Zero): Raw = 0.47947 -> 0.48
        expected_rti = [0.54 0.9 0.1 0.48]
        actual_rti = reef_tourism_index(rc, ce, sv, juv)

        @test actual_rti ≈ expected_rti
    end
end
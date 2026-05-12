using Test
using ADRIAIndicators

@testset "Coral Area Metrics" begin
    # 4D: [timesteps, groups, sizes, locations]
    # 5D: [timesteps, groups, sizes, locations, scenarios]
    
    n_timesteps = 2
    n_groups = 6
    n_sizes = 2
    n_locations = 3
    n_scenarios = 2
    
    # Use 0.05. Total cover per loc/step = 0.05 * 6 * 2 = 0.6.
    rel_cover_4d = fill(0.05, n_timesteps, n_groups, n_sizes, n_locations)
    habitat_area = [10.0, 20.0, 30.0] # km^2
    
    @testset "Absolute Coral Area 4D" begin
        abs_area = ADRIAIndicators.absolute_coral_area(rel_cover_4d, habitat_area)
        
        # Expected: 0.6 * habitat_area
        expected = 0.6 .* habitat_area'
        # Result should be (2, 3)
        @test size(abs_area) == (n_timesteps, n_locations)
        for t in 1:n_timesteps
            @test abs_area[t, :] ≈ expected[:]
        end
    end
    
    @testset "Absolute Coral Area 5D" begin
        rel_cover_5d = fill(0.05, n_timesteps, n_groups, n_sizes, n_locations, n_scenarios)
        # Scenario 2 has double cover
        rel_cover_5d[:, :, :, :, 2] .= 0.1
        
        abs_area = ADRIAIndicators.absolute_coral_area(rel_cover_5d, habitat_area)
        @test size(abs_area) == (n_timesteps, n_locations, n_scenarios)
        
        # Scenario 1: 0.6 * area
        # Scenario 2: 1.2 * area
        for t in 1:n_timesteps
            @test abs_area[t, :, 1] ≈ 0.6 .* habitat_area
            @test abs_area[t, :, 2] ≈ 1.2 .* habitat_area
        end
    end
    
    @testset "Coral Area Saved 4D vs 4D" begin
        # Intervention: 4D
        # Counterfactual: 4D
        intervention_cover = fill(0.1, n_timesteps, n_groups, n_sizes, n_locations)
        counterfactual_cover = fill(0.05, n_timesteps, n_groups, n_sizes, n_locations)
        
        saved_area = ADRIAIndicators.coral_area_saved(intervention_cover, counterfactual_cover, habitat_area)
        
        # Intervention Area = 1.2 * area
        # Counterfactual Area = 0.6 * area
        # Saved = 0.6 * area
        expected_saved = 0.6 .* habitat_area'
        
        @test size(saved_area) == (n_timesteps, n_locations)
        for t in 1:n_timesteps
            @test saved_area[t, :] ≈ expected_saved[:]
        end
    end
    
    @testset "Coral Area Saved 5D vs 5D" begin
        intervention_cover = fill(0.1, n_timesteps, n_groups, n_sizes, n_locations, n_scenarios)
        counterfactual_cover = fill(0.05, n_timesteps, n_groups, n_sizes, n_locations, n_scenarios)
        # Scenario 2 counterfactual has higher cover
        counterfactual_cover[:, :, :, :, 2] .= 0.08
        
        saved_area = ADRIAIndicators.coral_area_saved(intervention_cover, counterfactual_cover, habitat_area)
        
        @test size(saved_area) == (n_timesteps, n_locations, n_scenarios)
        
        # Scenario 1 Saved: (1.2 - 0.6) * area = 0.6 * area
        # Scenario 2 Saved: (1.2 - 0.96) * area = 0.24 * area
        # 0.08 * 6 * 2 = 0.96
        for t in 1:n_timesteps
            @test saved_area[t, :, 1] ≈ 0.6 .* habitat_area
            @test saved_area[t, :, 2] ≈ 0.24 .* habitat_area
        end
    end
end

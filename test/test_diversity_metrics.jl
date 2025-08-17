"""
Tests for coral_diversity and coral_evenness.
"""

using Test

using ReefMetrics: coral_diversity, coral_evenness

@testset "Coral Diversity" begin
    @testset "All Zero" begin
        rel_cover::Array{Float64,3} = zeros(Float64, 10, 5, 20)
        cor_div::Array{Float64,2} = coral_diversity(rel_cover)

        @test all(cor_div .== 0.0) || "A Matrix of all zeros should return 0 diversity."
    end

    @testset "One timestep, One Location" begin
        rel_cover::Array{Float64,3} = zeros(Float64, 1, 5, 1)
        rel_cover .= 0.2

        cor_div::Array{Float64,2} = coral_diversity(rel_cover)
        @test all(cor_div .≈ 0.8)

        rel_cover .= 0.1
        cor_div = coral_diversity(rel_cover)
        @test size(cor_div) == (1, 1)
        @test all(cor_div .≈ 0.8)

        rel_cover .= 0.05
        cor_div = coral_diversity(rel_cover)
        @test size(cor_div) == (1, 1)
        @test all(cor_div .≈ 0.8)

        rel_cover .= 0.0
        rel_cover[1, 1, 1] = 0.9
        cor_div = coral_diversity(rel_cover)
        @test size(cor_div) == (1, 1)
        @test all(cor_div .≈ 0.0)

        rel_cover .= 0.0
        rel_cover[1, 1, 1] = 0.8
        cor_div = coral_diversity(rel_cover)
        @test size(cor_div) == (1, 1)
        @test all(cor_div .≈ 0.0)

        rel_cover .= 0.0
        rel_cover[1, 1, 1] = 0.7
        cor_div = coral_diversity(rel_cover)
        @test size(cor_div) == (1, 1)
        @test all(cor_div .≈ 0.0)

        rel_cover .= 0.0
        rel_cover[1, :, 1] .= 0.0:0.1:0.4
        cor_div = coral_diversity(rel_cover)
        @test size(cor_div) == (1, 1)
        @test all(cor_div .≈ 0.7)
    end

    @testset "Multi-timestep, Multi-location" begin
        rel_cover::Array{Float64,3} = zeros(Float64, 3, 4, 3)
        rel_cover[1, :, 1] .= 0.25
        rel_cover[1, :, 2] .= 0.2
        rel_cover[1, :, 3] .= 0.15
        rel_cover[2, :, 1] .= 0.1:0.1:0.4
        rel_cover[2, :, 2] .= 0.0
        rel_cover[2, :, 3] .= 0.0:0.1:0.3
        rel_cover[3, 1, 1] = 0.1
        rel_cover[3, 2, 2] = 0.8
        rel_cover[3, 3, 3] = 0.3

        cor_div::Matrix{Float64} = coral_diversity(rel_cover)
        @test all(cor_div[1, 1] .≈ 0.75)
        @test all(cor_div[1, 2] .≈ 0.75)
        @test all(cor_div[1, 3] .≈ 0.75)
        @test all(cor_div[2, 1] .≈ 0.7)
        @test all(cor_div[2, 2] .≈ 0.0)
        @test all(cor_div[2, 3] .≈ 22 / 36)
        @test all(cor_div[3, 1] .≈ 0.0)
        @test all(cor_div[3, 2] .≈ 0.0)
        @test all(cor_div[3, 3] .≈ 0.0)
    end
end

@testset "Coral Evenness" begin
    @testset "All Zero" begin
        rel_cover = zeros(Float64, 10, 5, 20)
        evenness = coral_evenness(rel_cover)
        @test all(evenness .== 0.0)
    end

    @testset "Multi-timestep, Multi-location" begin
        n_tsteps, n_groups, n_locs = 2, 4, 2
        rel_cover = zeros(Float64, n_tsteps, n_groups, n_locs)

        rel_cover[1, :, 1] .= 0.2
        rel_cover[1, :, 2] = [0.4, 0.1, 0.1, 0.2]
        rel_cover[2, 1, 1] = 0.5

        evenness = coral_evenness(rel_cover)

        @test evenness[1, 1] ≈ 4.0
        @test evenness[1, 2] ≈ (1 / 0.34375)
        @test evenness[2, 1] ≈ 1.0
        @test evenness[2, 2] == 0.0

        out_evenness = coral_evenness(rel_cover)
        @test out_evenness == evenness
    end
end

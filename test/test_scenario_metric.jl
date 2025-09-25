using Test
using ADRIAIndicators: scenario_metric

@testset "Scenario Metric" begin
    metric_data = cat([1.0 2.0; 3.0 4.0], [5.0 6.0; 7.0 8.0], [9.0 10.0; 11.0 12.0], dims=3)
    # size: (2, 2, 3)

    location_area = [10.0, 20.0, 30.0]
    total_area = sum(location_area)
    location_dim = 3

    @testset "Relative In, Relative Out" begin
        result = scenario_metric(metric_data, location_area, location_dim)

        # Manual calculation
        abs_metric = metric_data .* reshape(location_area, (1, 1, 3))
        agg_metric = sum(abs_metric, dims=3)
        expected = agg_metric / total_area

        @test result ≈ dropdims(expected, dims=3)
    end

    @testset "Relative In, Absolute Out" begin
        result = scenario_metric(metric_data, location_area, location_dim; return_relative=false)

        # Manual calculation
        abs_metric = metric_data .* reshape(location_area, (1, 1, 3))
        expected = sum(abs_metric, dims=3)

        @test result ≈ dropdims(expected, dims=3)
    end

    @testset "Absolute In, Relative Out" begin
        result = scenario_metric(metric_data, location_area, location_dim; is_relative=false)

        # Manual calculation
        agg_metric = sum(metric_data, dims=3)
        expected = agg_metric / total_area

        @test result ≈ dropdims(expected, dims=3)
    end

    @testset "Absolute In, Absolute Out" begin
        result = scenario_metric(metric_data, location_area, location_dim; is_relative=false, return_relative=false)

        # Manual calculation
        expected = sum(metric_data, dims=3)

        @test result ≈ dropdims(expected, dims=3)
    end
end

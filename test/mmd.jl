using Test, Statistics, Distributions, Distances, DensityRatioEstimation
using DensityRatioEstimation: gaussian_gram, gaussian_gram_by_pairwise_sqd
using JuMP, Ipopt

@testset "MMDNumerical`" begin
    @testset "Positivity and normalisation constraints" begin
        dist_nu = Normal(1, 1)
        dist_de = Normal(0, 2)

        n_nu = 100
        x_nu = rand(dist_nu, 1, n_nu)
        n_de = 200
        x_de = rand(dist_de, 1, n_de)

        r_numerical = estimate_ratio(MMDNumerical(), x_de, x_nu)

        @test mean(r_numerical) ≈ 1
        @test all(r_numerical .> 0)
    end
end

@testset "MMDAnalytical`" begin
    @testset "Consistency between solve and inv methods" begin
        dist_nu = Normal(1, 1)
        dist_de = Normal(0, 2)

        n_nu = 100
        x_nu = rand(dist_nu, 1, n_nu)
        n_de = 200
        x_de = rand(dist_de, 1, n_de)

        r_solve = estimate_ratio(MMDAnalytical(method=:solve), x_de, x_nu)
        r_inv = estimate_ratio(MMDAnalytical(method=:inv), x_de, x_nu)

        @test r_solve ≈ r_inv

        r_solve = estimate_ratio(MMDAnalytical(0), x_de, x_nu)
        r_inv = estimate_ratio(MMDAnalytical(0), x_de, x_nu)

        @test r_solve ≈ r_inv
    end
end

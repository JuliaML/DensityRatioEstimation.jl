using Test, Statistics, Distributions, Distances, DensityRatioEstimation
using DensityRatioEstimation: pairwise_sqd, gaussian_gram, gaussian_gram_by_pairwise_sqd

@testset "Correctness of `pairwise_sqd` and `gaussian_gram_by_pairwise_sqd`" begin
    pairwise_sqd_golden(x) = pairwise(SqEuclidean(), x; dims=2)
    pairwise_sqd_golden(x, y) = pairwise(SqEuclidean(), x, y; dims=2)

    xtest = [
        1.0 2.0 4.0; 
        1.0 2.0 4.0
    ]

    ytest = [
        0.0 2.0 18.0;
        2.0 0.0  8.0;
        18.0 8.0  0.0
    ]

    @testset "`pairwise_sqd(x)`" begin
        n_randtests = 10
        @test pairwise_sqd(xtest) == ytest
        for _ = 1:n_randtests
            xrand = randn(784, 100)
            @test pairwise_sqd(xrand) ≈ pairwise_sqd_golden(xrand)
        end
    end

    @testset "`pairwise_sqd(x, y)`" begin
        n_randtests = 10
        @test pairwise_sqd(xtest, xtest) == ytest
        for _ = 1:n_randtests
            xrand = randn(784, 100)
            yrand = randn(784, 200)
            @test pairwise_sqd(xrand, yrand) ≈ pairwise_sqd_golden(xrand, yrand)
        end
    end

    @testset "`gaussian_gram`" begin
        for σ in sqrt.([1, 2, 4, 8, 16])
            K1 = gaussian_gram(xtest, σ)
            K2 = gaussian_gram(xtest, xtest, σ)
            K3 = gaussian_gram_by_pairwise_sqd(ytest, σ)
            @test K1 == K2 == K3 == exp.(-ytest / (2 * σ ^ 2))
        end
    end
end

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

        r_solve = estimate_ratio(MMDAnalytical(ϵ=0), x_de, x_nu)
        r_inv = estimate_ratio(MMDAnalytical(ϵ=0), x_de, x_nu)

        @test r_solve ≈ r_inv
    end
end
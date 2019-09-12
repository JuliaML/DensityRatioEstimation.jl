using Test, Statistics, Distributions, Distances, DensityRatioEstimation

@testset "Correctness of `pairwise_dot`" begin
    pairwise_dot(x) = pairwise(SqEuclidean(), x; dims=2)
    pairwise_dot(x, y) = pairwise(SqEuclidean(), x, y; dims=2)

    xtest = [
        1.0 2.0 4.0; 
        1.0 2.0 4.0
    ]

    ytest = [
        0.0 2.0 18.0;
        2.0 0.0  8.0;
        18.0 8.0  0.0
    ]

    @testset "`pairwise_dot(x)`" begin
        n_randtests = 10
        @test DensityRatioEstimation.pairwise_dot(xtest) == ytest
        for _ = 1:n_randtests
            xrand = randn(784, 100)
            @test DensityRatioEstimation.pairwise_dot(xrand) ≈ pairwise_dot(xrand)
        end
    end

    @testset "`pairwise_dot(x, y)`" begin
        n_randtests = 10
        @test DensityRatioEstimation.pairwise_dot(xtest) == ytest
        for _ = 1:n_randtests
            xrand = randn(784, 100)
            yrand = randn(784, 200)
            @test DensityRatioEstimation.pairwise_dot(xrand, yrand) ≈ pairwise_dot(xrand, yrand)
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

        r_numerical  = estimate_ratio(MMDNumerical(), x_de, x_nu)

        @test mean(r_numerical) ≈ 1
        @test all(r_numerical .> 0)
    end
end
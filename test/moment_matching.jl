using Test, Distances, DensityRatioEstimation

@testset "Correctness of `DensityRatioEstimation.pairwise_dot`" begin
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

    @testset "`DensityRatioEstimation.pairwise_dot(x)`" begin
        n_randtests = 10
        @test DensityRatioEstimation.pairwise_dot(xtest) == ytest
        for _ = 1:n_randtests
            xrand = randn(784, 100)
            @test DensityRatioEstimation.pairwise_dot(xrand) ≈ pairwise_dot(xrand)
        end
    end

    @testset "`DensityRatioEstimation.pairwise_dot(x, y)`" begin
        n_randtests = 10
        @test DensityRatioEstimation.pairwise_dot(xtest) == ytest
        for _ = 1:n_randtests
            xrand = randn(784, 100)
            yrand = randn(784, 200)
            @test DensityRatioEstimation.pairwise_dot(xrand, yrand) ≈ pairwise_dot(xrand, yrand)
        end
    end
end

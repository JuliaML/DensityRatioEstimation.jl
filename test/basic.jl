@testset "Basic" begin
  @testset "Gramian" begin
    x_nu, x_de = [rand(2) for i=1:100], [rand(2) for i=1:200]
    G = DensityRatioEstimation.gaussian_gramian(x_nu, x_de, σ=1.0)
    @test size(G) == (length(x_nu), length(x_de))
    @test all(G .> 0)

    G = DensityRatioEstimation.gaussian_gramian(x_nu, x_nu, σ=2.0)
    @test issymmetric(G)
    @test all(G .> 0)
  end
end

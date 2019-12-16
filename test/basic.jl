@testset "Basic" begin
  for (d_nu, d_de) in [pair₁, pair₂]
    Random.seed!(123)
    x_nu, x_de = rand(d_nu, 100), rand(d_de, 200)

    @testset "Gramian" begin
      G = DensityRatioEstimation.gaussian_gramian(x_nu, x_de, σ=1.0)
      @test size(G) == (length(x_nu), length(x_de))
      @test all(G .> 0)

      G = DensityRatioEstimation.gaussian_gramian(x_nu, x_nu, σ=2.0)
      @test issymmetric(G)
      @test all(G .> 0)
    end

    @testset "$dre -- $optlib" for (dre, optlib) in [(KMM(),   JuMPLib),
                                                     (KLIEP(), OptimLib),
                                                     (KLIEP(), ConvexLib)]

      r = densratio(x_nu, x_de, dre, optlib=optlib)

      # density ratios must be positive
      @test all(r .> 0)
    end
  end
end

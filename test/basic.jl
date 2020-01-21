@testset "Basic" begin
  @testset "Gramian" begin
    x_nu, x_de = [rand(2) for i=1:100], [rand(2) for i=1:200]
    G = DensityRatioEstimation.gaussian_gramian(x_nu, x_de, σ=1.0)
    @test size(G) == (length(x_nu), length(x_de))
    @test all(G .> 0)

    G = DensityRatioEstimation.gaussian_gramian(x_nu, x_nu, σ=2.0)
    @test issymmetric(G)
    @test all(G .> 0)

    # features can be any indexable
    x_nu = [(a=1.,b=2.),(a=3.,b=4.)]
    x_de = [(a=1.,b=2.),(a=3.,b=4.),(a=5.,b=6.)]
    G = DensityRatioEstimation.gaussian_gramian(x_nu, x_de)
    @test size(G) == (2, 3)
    @test all(G .> 0)
  end

  for (d_nu, d_de) in [pair₁, pair₂]
    Random.seed!(123)
    x_nu, x_de = rand(d_nu, 100), rand(d_de, 200)
    @testset "$dre -- $optlib" for (dre, optlib) in [(KMM(),   JuMPLib),
                                                     (KLIEP(), OptimLib),
                                                     (KLIEP(), ConvexLib),
                                                     (LSIF(), JuMPLib)]

      r = densratio(x_nu, x_de, dre, optlib=optlib)

      # density ratios must be positive
      @test all(r .> 0)
    end
  end
end

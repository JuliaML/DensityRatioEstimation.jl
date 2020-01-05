@testset "Basic" begin 
  for (i, pair) in enumerate([pair₁, pair₂])
    d_nu, d_de = pair
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

    @testset "$drestr -- $optlib" for (drestr, dre, optlib) in [
      ("KMM",	  KMM(σ=1.5, B=Inf, ϵ=0.01, λ=0.001), JuMPLib),
      ("KLIEP", KLIEP(), ConvexLib)
    ]
      r̂ = densratio(x_nu, x_de, dre, optlib=optlib)

      # density ratios must be positive
      @test all(r̂ .> 0)

      if i == 1
        # compare against true ratio
        r = pdf.(d_nu, x_de) ./ pdf.(d_de, x_de)
        @test r ≈ r̂ rtol=2e-1
      end
    end
  end
end

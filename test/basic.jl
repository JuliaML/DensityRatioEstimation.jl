@testset "Basic" begin
  for (d_nu, d_de) in [pair₁, pair₂]
    Random.seed!(123)
    x_nu, x_de = rand(d_nu, 100), rand(d_de, 200)
    @testset "$dre -- $optlib" for (dre, optlib) in [(KMM(),   JuMPLib),
                                                     (KLIEP(), OptimLib),
                                                     (KLIEP(), ConvexLib)]

      r = densratio(x_nu, x_de, dre, optlib=optlib)

      # density ratios must be positive
      @test all(r .> 0)
    end
  end
end

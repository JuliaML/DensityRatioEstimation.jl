@testset "KLIEP -- $optlib" for optlib in [OptimLib, ConvexLib]
  for (pair, rtol) in [(pair₁, 2e-1), (pair₂, 4e-1)]
    d_nu, d_de = pair
    Random.seed!(123)
    x_nu, x_de = rand(d_nu, 1_000), rand(d_de, 500)

    # estimated density ratio
    σ, b = 1.0, 100
    r̂ = densratio(x_nu, x_de, KLIEP(σ=σ, b=b), optlib=optlib)

    # simplex constraints
    @test abs(mean(r̂) - 1) ≤ 1e-2
    @test all(r̂ .≤ Inf)

    r = pdf.(d_nu, x_de) ./ pdf.(d_de, x_de)
    @test r ≈ r̂ rtol=rtol

    if visualtests
      gr(size=(800, 800))
      @plottest plot_d_nu(pair, x_de, r̂) joinpath(datadir, "KLIEP-$optlib-$i.png") !istravis
    end
  end
end

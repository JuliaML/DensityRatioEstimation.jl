@testset "KLIEP -- $optlib" for optlib in [OptimLib, ConvexLib]
  for (i, (pair, rtol)) in enumerate([(pair₁, 5e-1), (pair₂, 5e-1)])
    d_nu, d_de = pair
    rng = StableRNG(42)
    x_nu, x_de = rand(rng, d_nu, 1_000), rand(rng, d_de, 500)

    # estimated density ratio
    kliep = KLIEP(σ=1.0, b=100, rng=rng)
    r̂ = densratio(x_nu, x_de, kliep, optlib=optlib)

    # simplex constraints
    @test abs(mean(r̂) - 1) ≤ 1e-2
    @test all(r̂ .≤ Inf)

    r = pdf.(d_nu, x_de) ./ pdf.(d_de, x_de)
    @test r ≈ r̂ rtol = rtol

    if visualtests
      @test_reference "data/KLIEP-$optlib-$i.png" plot_d_nu(pair, x_de, r̂)
    end
  end
end

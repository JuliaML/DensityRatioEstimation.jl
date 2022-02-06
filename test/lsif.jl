@testset "LSIF -- $optlib" for optlib in [OptimLib, JuMPLib]
  for (i, (pair, rtol)) in enumerate([(pair₁, 2e-1), (pair₂, 4e-1)])
    d_nu, d_de = pair
    rng = MersenneTwister(42)
    x_nu, x_de = rand(rng, d_nu, 1000), rand(rng, d_de, 500)

    # estimated density ratio
    lsif = LSIF(σ=1.0, b=100)
    r̂ = densratio(x_nu, x_de, lsif, optlib=optlib)

    # simplex constraints
    @test abs(mean(r̂) - 1) ≤ 2e-2
    @test all(r̂ .≤ Inf)

    r = pdf.(d_nu, x_de) ./ pdf.(d_de, x_de)
    @test r ≈ r̂ rtol=rtol

    if visualtests
      gr(size=(800,800))
      @test_reference "data/LSIF-$optlib-$i.png" plot_d_nu(pair,x_de,r̂) 
    end
  end
end

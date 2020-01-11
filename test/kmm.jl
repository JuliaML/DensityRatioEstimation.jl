@testset "KMM -- $optlib" for optlib in [JuliaLib, JuMPLib]
  for (i, pair) in enumerate([pair₁, pair₂])
    d_nu, d_de = pair
    Random.seed!(123)
    x_nu, x_de = rand(d_nu, 2_000), rand(d_de, 1_000)

    # estimated density ratio
    σ, B, ϵ, λ = 1.5, Inf, 0.01, 0.001
    r̂ = densratio(x_nu, x_de, KMM(σ=σ, B=B, ϵ=ϵ,  λ=λ), optlib=optlib)

    # simplex constraints
    @test abs(mean(r̂) - 1) ≤ 1e-2
    @test all(r̂ .≤ B)

    if i == 1 # FIXME: only check correctness for the Gaussian case now
      # compare against true ratio
      r = pdf.(d_nu, x_de) ./ pdf.(d_de, x_de)
      @test r ≈ r̂ rtol=2e-1
    end

    if visualtests
      gr(size=(800, 800))
      @plottest plot_d_nu(pair, x_de, r̂) joinpath(datadir, "KMM-$optlib-$i.png") !istravis
    end
  end
end

@testset "$(nameof(dreType)) -- $optlib" for (dreType, optlib) in zip([uKMM, uKMM, KMM], [JuliaLib, JuMPLib, JuMPLib])
  for (i, (pair, rtol)) in enumerate([(pair₁, 2e-1), (pair₂, 4e-1)])
    d_nu, d_de = pair
    rng = MersenneTwister(42)
    x_nu, x_de = rand(rng, d_nu, 2_000), rand(rng, d_de, 1_000)

    # estimated density ratio
    D = [sqrt(DensityRatioEstimation.euclidsq(x, y)) for x in x_nu, y in x_de]
    σ, λ = median(D), 0.01
    kmm = dreType(σ=σ, λ=λ)
    r̂ = densratio(x_nu, x_de, kmm; optlib=optlib)

    @test abs(mean(r̂) - 1) ≤ 1e-2

    # simplex constraints
    if dreType == KMM
      @test all(0 .≤ r̂)
      @test all(r̂ .≤ kmm.B)
      @test all((1-kmm.ϵ) ≤ mean(r̂) ≤ (1+kmm.ϵ))
    end

    # compare against true ratio
    r = pdf.(d_nu, x_de) ./ pdf.(d_de, x_de)
    @test r ≈ r̂ rtol=rtol

    # type consistency
    @test eltype(r) == typeof(σ)

    # iterator and matrix version consistency for JuliaLib
    if optlib == JuliaLib
      r̂_mat = densratio(reshape(x_nu, 1, :), reshape(x_de, 1, :), kmm; optlib=optlib)
      @test r̂ ≈ r̂_mat
    end

    if visualtests
      gr(size=(800, 800))
      @test_reference "data/$(nameof(dreType))-$optlib-$i.png" plot_d_nu(pair, x_de, r̂)
    end
  end
end

@testset "KMM -- $optlib" for optlib in [JuliaLib, JuMPLib]
  for (i, (pair, rtol)) in enumerate([(pair₁, 2e-1), (pair₂, 4e-1)])
    d_nu, d_de = pair
    rng = MersenneTwister(42)
    x_nu, x_de = rand(rng, d_nu, 2_000), rand(rng, d_de, 1_000)

    # estimated density ratio
    D = [sqrt(DensityRatioEstimation.euclidsq(x, y)) for x in x_nu, y in x_de]
    σ, B, ϵ, λ = median(D), Inf, 0.001, 0.01
    kmm = KMM(σ=σ, B=B, ϵ=ϵ,  λ=λ)
    r̂ = densratio(x_nu, x_de, kmm; optlib=optlib)

    # simplex constraints
    @test abs(mean(r̂) - 1) ≤ 1e-2
    @test all(r̂ .≤ B)

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
      @test_reference "data/KMM-$optlib-$i.png" plot_d_nu(pair, x_de, r̂)
    end
  end
end

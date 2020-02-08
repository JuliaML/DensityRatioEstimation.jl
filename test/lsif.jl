@testset "LSIF" begin
  for (i, (pair, rtol)) in enumerate([(pair₁, 2e-1), (pair₂, 4e-1)])
    d_nu, d_de = pair
    Random.seed!(123)
    x_nu, x_de = rand(d_nu, 1000), rand(d_de, 500)

    # estimated density ratio
    r̂ = densratio(x_nu, x_de, LSIF(σ=1.0, b=100))

    # simplex constraints
    @test abs(mean(r̂) - 1) ≤ 2e-2
    @test all(r̂ .≤ Inf)

    r = pdf.(d_nu, x_de) ./ pdf.(d_de, x_de)
    @test r ≈ r̂ rtol=rtol

    if visualtests
      gr(size=(800,800))
      @plottest plot_d_nu(pair,x_de,r̂) joinpath(datadir,"LSIF-$i.png") !istravis
    end
  end
end

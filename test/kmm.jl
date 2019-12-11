@testset "KMM" begin
  d_nu, d_de = pair₁

  Random.seed!(123)
  x_nu, x_de = rand(d_nu, 2000), rand(d_de, 1000)

  # estimated density ratio
  σ, B, ϵ = 1.0, 100.0, 0.01
  r̂ = densratio(x_nu, x_de, KMM(σ, B, ϵ))

  # simplex constraints
  @test abs(mean(r̂) - 1) ≤ ϵ
  @test all(r̂ .≤ B)

  if visualtests
    gr(size=(800,800))
    @plottest plot_d_nu(pair₁,x_de,r̂) joinpath(datadir,"KMM.png") !istravis
  end
end

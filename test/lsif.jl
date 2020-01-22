@testset "LSIF" begin
  for (i, pair) in enumerate([pair₁, pair₂])
    d_nu, d_de = pair
    Random.seed!(123)
    x_nu, x_de = rand(d_nu, 1000), rand(d_de, 500)

    # estimated density ratio
    r̂ = densratio(x_nu, x_de, LSIF(σ=1.0, b=100))

    if visualtests
      gr(size=(800,800))
      @plottest plot_d_nu(pair,x_de,r̂) joinpath(datadir,"LSIF-$i.png") !istravis
    end
  end
end

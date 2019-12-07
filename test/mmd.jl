@testset "MMD" begin
  n_nu, n_de = 100, 200
  dist_nu = Normal(1, 1)
  dist_de = Normal(0, 2)
  x_nu = rand(dist_nu, 1, n_nu)
  x_de = rand(dist_de, 1, n_de)

  @testset "Julia" begin
    r_solve = estimate_ratio(MMDAnalytical(method=:solve), x_de, x_nu)
    r_inv = estimate_ratio(MMDAnalytical(method=:inv), x_de, x_nu)

    @test r_solve ≈ r_inv

    r_solve = estimate_ratio(MMDAnalytical(0), x_de, x_nu)
    r_inv = estimate_ratio(MMDAnalytical(0), x_de, x_nu)

    @test r_solve ≈ r_inv
  end

  @testset "JuMP" begin
    r = estimate_ratio(MMDNumerical(), x_de, x_nu)
    @test mean(r) ≈ 1
    @test all(r .> 0)
  end
end

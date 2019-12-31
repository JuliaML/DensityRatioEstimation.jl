@testset "KLIEP" begin
	for (i, pair) in enumerate([pair₁, pair₂])
		d_nu, d_de = pair
		Random.seed!(123)
		x_nu, x_de = rand(d_nu, 1_000), rand(d_de, 500)

		# estimated density ratio
		σ, b = 1.0, 100
		r̂ = densratio(x_nu, x_de, KLIEP(σ=σ, b=b))

		if i == 1
			# compare against true ratio
			r = pdf.(d_nu, x_de) ./ pdf.(d_de, x_de)
			@test r ≈ r̂ rtol=2e-1
		end

		if visualtests
			gr(size=(800, 800))
			@plottest plot_d_nu(pair, x_de, r̂) joinpath(datadir, "KLIEP-$i.png") !istravis
		end
	end
end

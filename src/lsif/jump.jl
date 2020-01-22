# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

using .JuMP
using .Ipopt

function _lsif_coeffs(x_nu, x_de, centers::AbstractVector{Int},
                      dre::LSIF, optlib::Type{JuMPLib})
  # retrieve parameters
  σ, λ, b = dre.σ, dre.λ, length(centers)

  # number of numerator and denominator samples
  n_nu, n_de = length(x_nu), length(x_de)

  # constants for objective and constraints
  K_nu = gaussian_gramian(x_nu, x_nu[centers], σ=σ)
  K_de = gaussian_gramian(x_de, x_nu[centers], σ=σ)
  Φ = Matrix{eltype(K_de)}(undef, b, b)
  for l′ in 1:b
    φ′ = view(K_de, :, l′)
    for l in 1:l′
      φ = view(K_de, :, l)
      Φ[l,l′] = mean(φ .* φ′)
    end
  end
  H = Symmetric(Φ)
  h = vec(mean(K_nu, dims=1))

  model = Model(with_optimizer(Ipopt.Optimizer, print_level=0, sb="yes"))
  @variable(model, α[1:b])
  @objective(model, Min, (1/2) * dot(α, H*α - 2h) + λ * sum(α))
  @constraint(model, α .≥ 0)

  # solve the problem
  optimize!(model)

  # density ratio
  value.(α)
end

# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

using .Optim

function _lsif_coeffs(x_nu, x_de, centers::AbstractVector{Int},
                      dre::LSIF, optlib::Type{OptimLib})
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

  # constants for inequality constraints
  T = eltype(x_de[1])
  lx = fill(zero(T), b)
  ux = fill(Inf, b)

  # objective
  f(α)       = (1/2) * dot(α, H*α - 2h) + λ * sum(α)
  ∇f!(grad, α)  = (grad .= H*α .- h .+ λ)
  ∇²f!(hess, α) = (hess .= H)

  # initial guess
  αₒ = fill(one(T), b)

  # optimization problem
  objective   = TwiceDifferentiable(f, ∇f!, ∇²f!, αₒ)
  constraints = TwiceDifferentiableConstraints(lx, ux)
  initguess   = αₒ

  # solve problem with interior-point primal-dual Newton
  solution = optimize(objective, constraints, initguess, IPNewton())

  # optimal coefficients
  solution.minimizer
end

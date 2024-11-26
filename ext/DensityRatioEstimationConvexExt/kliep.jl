# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function DensityRatioEstimation._kliep_coeffs(K_nu, K_de, dre::KLIEP, optlib::Type{ConvexLib})
  # retrieve parameters
  σ, b = dre.σ, size(K_de, 2)

  # number of numerator and denominator samples
  n_nu, n_de = size(K_nu, 1), size(K_de, 1)

  # constants for objective and constraints
  K = K_nu
  k = sum(K_de, dims=1)

  # objective function and constraints
  α = Convex.Variable(b)
  w = K * α
  objective = sum(log(w[i]) for i in 1:n_nu)
  constraints = [α ≥ 0, dot(α, k) == n_de]

  # optimization problem
  problem = Convex.maximize(objective, constraints)

  # solve problem with ECOS solver
  Convex.solve!(problem, ECOS.Optimizer, silent=true)

  # optimal coefficients
  vec(α.value)
end

# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

import .Convex
import .ECOS: ECOSSolver

function _densratio(x_nu, x_de, dre::KLIEP, optlib::Type{ConvexLib})
  # retrieve parameters
  @unpack σ, b = dre

  # number of numerator and denominator samples
  n_nu, n_de = length(x_nu), length(x_de)

  @assert b ≤ n_nu "more basis elements than numerator samples"

  # basis for kernel approximation
  basis = sample(1:n_nu, b, replace=false)

  # constants for objective and constraints
  K = gaussian_gramian(x_nu, x_nu[basis], σ=σ)
  P = gaussian_gramian(x_de, x_nu[basis], σ=σ)
  p = sum(P, dims=1)

  # objective function and constraints
  α = Convex.Variable(b); w = K*α
  objective   = sum(log(w[i]) for i in 1:n_nu)
  constraints = [α ≥ 0, dot(α,p) == n_de]

  # optimization problem
  problem = Convex.maximize(objective, constraints)

  # solve problem with ECOS solver
  Convex.solve!(problem, ECOSSolver(verbose=false))

  # optimal coefficients
  a = vec(α.value)

  # density ratio
  P*a
end

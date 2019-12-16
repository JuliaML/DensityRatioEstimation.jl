# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

import .Convex
import .ECOS: ECOSSolver

function _kliep_coeffs(x_nu, x_de, centers::AbstractVector{Int},
                       dre::KLIEP, optlib::Type{ConvexLib})
  # retrieve parameters
  @unpack σ, b = dre

  # number of numerator and denominator samples
  n_nu, n_de = length(x_nu), length(x_de)

  # constants for objective and constraints
  K = gaussian_gramian(x_nu, x_nu[centers], σ=σ)
  P = gaussian_gramian(x_de, x_nu[centers], σ=σ)
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
  vec(α.value)
end

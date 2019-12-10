# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

import .Convex
import .ECOS: ECOSSolver

function _density_ratio(x_nu, x_de, dre::KLIEP, optlib::Type{ConvexLib})
  # retrieve parameters
  σ, b = dre.σ, dre.b

  # number of numerator and denominator samples
  n_nu = length(x_nu)
  n_de = length(x_de)

  # basis for kernel approximation
  basis = sample(1:n_nu, b, replace=false)

  # constants for objective and constraints
  K = gaussian_gramian(x_nu, x_nu[basis], σ=σ)
  P = gaussian_gramian(x_de, x_nu[basis], σ=σ)
  p = sum(P, dims=1)

  # objective function
  α = Convex.Variable(b); w = K*α
  J = sum(log(w[i]) for i in 1:n_nu)

  # optimization problem
  problem = Convex.maximize(J, [α ≥ 0, dot(α,p) == n_de])

  # solve problem with ECOS solver
  Convex.solve!(problem, ECOSSolver(verbose=false))

  # density ratio
  r = P*vec(α.value)
end

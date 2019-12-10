# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

using .JuMP
using .Ipopt

function _density_ratio(x_nu, x_de, dre::KMM, optlib::Type{JuMPLib})
  # retrieve parameters
  σ = dre.σ

  # number of numerator and denominator samples
  n_nu = length(x_nu)
  n_de = length(x_de)

  # constants for objective and constraints
  Kdede = gaussian_gramian(x_de, x_de, σ=σ)
  Kdenu = gaussian_gramian(x_de, x_nu, σ=σ)

  # optimization problem
  model = Model(with_optimizer(Ipopt.Optimizer, print_level=0, sb="yes"))
  @variable(model, r[1:n_de])
  @objective(model, Min, 1 / n_de ^ 2 * sum(r[i] * Kdede[i,j] * r[j] for i = 1:n_de, j=1:n_de) - 2 / (n_de * n_nu) * sum(r[i] * Kdenu[i,j] for i = 1:n_de, j=1:n_nu))
  @constraint(model, r .>= 0)
  @constraint(model, 1 / n_de * sum(r) == 1)

  # solve the optimization problem
  optimize!(model)

  # density ratio
  value.(r)
end

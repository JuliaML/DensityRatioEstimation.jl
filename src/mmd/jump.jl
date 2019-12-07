# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

using .JuMP
using .Ipopt

struct MMDNumerical <: AbstractMMD
  positivity::Bool
  normalisation::Bool
end

MMDNumerical(;positivity::Bool=true, normalisation::Bool=true) =
  MMDNumerical(positivity, normalisation)

function _density_ratio(mmd::MMDNumerical, Kdede, Kdenu)
  n_de, n_nu = size(Kdenu)

  model = Model(with_optimizer(Ipopt.Optimizer, print_level=0, sb="yes"))
  @variable(model, r[1:n_de])
  @objective(model, Min, 1 / n_de ^ 2 * sum(r[i] * Kdede[i,j] * r[j] for i = 1:n_de, j=1:n_de) - 2 / (n_de * n_nu) * sum(r[i] * Kdenu[i,j] for i = 1:n_de, j=1:n_nu))
  if mmd.positivity
    @constraint(model, r .>= 0)
  end
  if mmd.normalisation
    @constraint(model, 1 / n_de * sum(r) == 1)
  end

  optimize!(model)

  value.(r)
end

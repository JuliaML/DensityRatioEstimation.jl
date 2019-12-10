# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

using .JuMP
using .Ipopt

function _density_ratio(x_nu, x_de, dre::KMM, optlib::Type{JuMPLib})
  # retrieve parameters
  σ, B = dre.σ, dre.B

  # number of numerator and denominator samples
  m′, m = length(x_nu), length(x_de)

  # ϵ ∼ O(B/√m) (Huang et al. 2006.)
  ϵ = B / √m

  # constants for objective and constraints
  K = gaussian_gramian(x_de, x_de, σ=σ)
  A = gaussian_gramian(x_de, x_nu, σ=σ)
  κ = (m / m′) * sum(A, dims=2)

  # optimization problem
  model = Model(with_optimizer(Ipopt.Optimizer, print_level=0, sb="yes"))
  @variable(model, β[1:m])
  @objective(model, Min, (1/2)*β'*K*β - κ'*β)
  @constraint(model, 0 .≤ β .≤ B)
  @constraint(model, (1-ϵ)m ≤ sum(β) ≤ (1+ϵ)m)

  # solve the optimization problem
  optimize!(model)

  # density ratio
  value.(β)
end

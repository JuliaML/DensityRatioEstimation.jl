# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENSE in the project root.
# ------------------------------------------------------------------

using .JuMP
using .Ipopt

function _kmm_ratios(K, κ, dre::T, optlib::Type{JuMPLib}) where {T<:AbstractFloat}
  # retrieve parameters
  @unpack B, ϵ = dre

  # optimization problem
  model = Model(with_optimizer(Ipopt.Optimizer, print_level=0, sb="yes"))
  @variable(model, β[1:m])
  @objective(model, Min, (1/2) * dot(β, K*β - 2κ))
  @constraint(model, 0 .≤ β)
  isinf(B) || @constraint(model, β .≤ B)
  @constraint(model, (1-ϵ)m ≤ sum(β) ≤ (1+ϵ)m)

  # solve the problem
  optimize!(model)

  # density ratio
  value.(β)
end

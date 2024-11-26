# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# This file is part of the module DensityRatioEstimationJuMPExt.

function _kmm_jump_model(K, κ, dre::AbstractKMM, optlib::Type{JuMPLib})
  # number of denominator samples
  m = length(κ)

  # optimization problem
  model = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" => "yes"))
  @variable(model, β[1:m])
  @objective(model, Min, (1 / 2) * dot(β, K * β - 2κ))

  return model, β
end

function DensityRatioEstimation._kmm_ratios(K, κ, dre::uKMM, optlib::Type{JuMPLib})
  # build the problem without constraints
  model, β = _kmm_jump_model(K, κ, dre, optlib)

  # solve the problem
  optimize!(model)

  # density ratio
  value.(β)
end

function DensityRatioEstimation._kmm_ratios(K, κ, dre::KMM, optlib::Type{JuMPLib})
  # retrieve parameters
  @unpack B, ϵ = dre

  # build the problem without constraints
  model, β = _kmm_jump_model(K, κ, dre, optlib)

  # adding constraints
  @constraint(model, 0 .≤ β)
  isinf(B) || @constraint(model, β .≤ B)
  @constraint(model, (1 - ϵ) ≤ mean(β) ≤ (1 + ϵ))

  # solve the problem
  optimize!(model)

  # density ratio
  value.(β)
end

# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# This file is part of the module DensityRatioEstimationJuMPExt.

function DensityRatioEstimation._lsif_coeffs(H, h, dre::LSIF, optlib::Type{JuMPLib})
  # retrieve parameters
  λ, b = dre.λ, length(h)

  model = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" => "yes"))
  @variable(model, α[1:b])
  @objective(model, Min, (1 / 2) * dot(α, H * α - 2h) + λ * sum(α))
  @constraint(model, α .≥ 0)

  # solve the problem
  optimize!(model)

  # density ratio
  value.(α)
end

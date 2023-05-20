# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module DensityRatioEstimationChainRulesCoreExt

if isdefined(Base, :get_extension)
  using DensityRatioEstimation
  using ChainRulesCore
else
  using ..DensityRatioEstimation
  using ..ChainRulesCore
end

ChainRulesCore.@non_differentiable DensityRatioEstimation.safe_diagm(::Any, ::Any)

end #module

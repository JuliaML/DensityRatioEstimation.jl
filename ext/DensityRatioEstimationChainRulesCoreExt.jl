# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module DensityRatioEstimationChainRulesCoreExt

using DensityRatioEstimation
using ChainRulesCore

ChainRulesCore.@non_differentiable DensityRatioEstimation.safe_diagm(::Any, ::Any)

end #module

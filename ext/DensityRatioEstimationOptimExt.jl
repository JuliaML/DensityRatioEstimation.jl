# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module DensityRatioEstimationOptimExt

if isdefined(Base, :get_extension)
  using DensityRatioEstimation
  using DensityRatioEstimation: KLIEP, LSIF, OptimLib
  using Optim
else
  using ..DensityRatioEstimation
  using ..DensityRatioEstimation: KLIEP, LSIF, OptimLib
  using ..Optim
end

include("../src/kliep/optim.jl")
include("../src/lsif/optim.jl")

end #module

# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module DensityRatioEstimationConvexExt

if isdefined(Base, :get_extension)
  using DensityRatioEstimation
  using DensityRatioEstimation: KLIEP, ConvexLib
  using Convex
  using ECOS
else
  using ..DensityRatioEstimation
  using ..DensityRatioEstimation: KLIEP, ConvexLib
  using ..Convex
  using ..ECOS
end

include("../src/kliep/convex.jl")

end #module

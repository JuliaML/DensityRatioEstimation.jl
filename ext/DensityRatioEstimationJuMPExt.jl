# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module DensityRatioEstimationJuMPExt

if isdefined(Base, :get_extension)
  using DensityRatioEstimation
  using DensityRatioEstimation: LSIF, JuMPLib, AbstractKMM, uKMM, KMM
  using DensityRatioEstimation.Parameters
  using JuMP
  using Ipopt
  using LinearAlgebra
  using Statistics
else
  using ..DensityRatioEstimation
  using ..DensityRatioEstimation: LSIF, JuMPLib, AbstractKMM, uKMM, KMM
  using ..DensityRatioEstimation.Parameters
  using ..JuMP
  using ..Ipopt
  using ..LinearAlgebra
  using ..Statistics
end

include("../src/kmm/jump.jl")
include("../src/lsif/jump.jl")

end #module

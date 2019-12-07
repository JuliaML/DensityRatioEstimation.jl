# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module DensityRatioEstimation

using Statistics
using LinearAlgebra
using Distances
using Parameters # TODO: eliminate dependency

using Requires
function __init__()
  # load methods based on available dependencies
  @require JuMP="4076af6c-e467-56ae-b986-b466b2749572" begin
    @require Ipopt="b6b21f68-93f8-5de0-b562-5493be1d77c9" include("mmd/jump.jl")
  end
end

include("mmd/julia.jl")

export
  # types
  DensityRatioEstimator,
  MMD,
  MMDAnalytical, # TODO: deprecate
  MMDNumerical, # TODO: deprecate

  # functions
  density_ratio

end # module

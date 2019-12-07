# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module DensityRatioEstimation

using Statistics
using LinearAlgebra
using Distances

# API for density ratio estimation
include("api.jl")

# pure Julia implementations
include("mmd/julia.jl")

# implementations with extra dependencies
using Requires
function __init__()
  @require JuMP="4076af6c-e467-56ae-b986-b466b2749572" begin
    @require Ipopt="b6b21f68-93f8-5de0-b562-5493be1d77c9" include("mmd/jump.jl")
  end
end

export
  # types
  DensityRatioEstimator,
  MMD,
  MMDAnalytical, # TODO: deprecate
  MMDNumerical, # TODO: deprecate

  # functions
  density_ratio

end # module

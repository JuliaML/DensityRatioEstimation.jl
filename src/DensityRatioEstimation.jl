# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module DensityRatioEstimation

using StatsBase
using Statistics
using LinearAlgebra
using Parameters

# implement fit for estimators
import StatsBase: fit

# API for density ratio estimation
include("api/optlibs.jl")
include("api/estimators.jl")

# utility functions
include("utils.jl")

# available estimators
include("kmm.jl")
include("kliep.jl")

# pure Julia implementations
include("kmm/julia.jl")

# implementations that require extra dependencies
using Requires
function __init__()
  # KMM
  @require JuMP="4076af6c-e467-56ae-b986-b466b2749572" begin
    @require Ipopt="b6b21f68-93f8-5de0-b562-5493be1d77c9" include("kmm/jump.jl")
  end

  # KLIEP
  @require Optim="429524aa-4258-5aef-a3af-852621145aeb" include("kliep/optim.jl")
  @require Convex="f65535da-76fb-5f13-bab9-19810c17039a" begin
    @require ECOS="e2685f51-7e38-5353-a97d-a921fd2c8199" include("kliep/convex.jl")
  end
end

export
  # estimators
  DensityRatioEstimator,
  KMM, KLIEP,

  # optim libs
  JuliaLib,
  OptimLib,
  ConvexLib,
  JuMPLib,

  # functions
  densratio,
  fit

end # module

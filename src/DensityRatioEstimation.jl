# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENSE in the project root.
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
include("api/fitters.jl")

# utility functions
include("utils.jl")

# available estimators
include("kmm.jl")
include("kliep.jl")
include("lsif.jl")

# available estimator fitters
include("lcv.jl")

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

  # LSIF
  @require Optim="429524aa-4258-5aef-a3af-852621145aeb" include("lsif/optim.jl")
  @require JuMP="4076af6c-e467-56ae-b986-b466b2749572" begin
    @require Ipopt="b6b21f68-93f8-5de0-b562-5493be1d77c9" include("lsif/jump.jl")
  end

  # AD and GPU libs
  @require Zygote="e88e6eb3-aa80-5325-afca-941959d7151f" include("lib/zygote.jl")
  @require CuArrays="3a865a2d-5b23-5a0f-bc46-62713ec82fae" include("lib/cuarrays.jl")
end

export
  # optim libs
  OptimizationLibrary,
  JuliaLib,
  OptimLib,
  ConvexLib,
  JuMPLib,

  # estimators
  DensityRatioEstimator,
  KMM, KLIEP, LSIF,
  available_optlib,
  default_optlib,
  densratiofunc,
  densratio,

  # estimator fitters
  EstimatorFitter,
  LCV,
  fit

end # module

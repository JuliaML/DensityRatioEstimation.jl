# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module DensityRatioEstimation

using StatsBase
using Statistics
using LinearAlgebra
using Parameters
using Random

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

if !isdefined(Base, :get_extension)
  using Requires
end
# implementations that require extra dependencies
@static if !isdefined(Base, :get_extension)
  function __init__()

    #Solvers

    #JuMP: KMM, LSIF
    @require JuMP = "4076af6c-e467-56ae-b986-b466b2749572" begin
      @require Ipopt = "b6b21f68-93f8-5de0-b562-5493be1d77c9" begin
        include("../ext/DensityRatioEstimationJuMPExt.jl")
      end
    end
    #Optim: KLIEP, LSIF
    @require Optim = "429524aa-4258-5aef-a3af-852621145aeb" begin
      include("../ext/DensityRatioEstimationOptimExt.jl")
    end

    #Convex: KLIEP
    @require Convex = "f65535da-76fb-5f13-bab9-19810c17039a" begin
      @require ECOS = "e2685f51-7e38-5353-a97d-a921fd2c8199" begin
        include("../ext/DensityRatioEstimationConvexExt.jl")
      end
    end

    # AD and GPU libs
    @require ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4" begin
      include("../ext/DensityRatioEstimationChainRulesCoreExt.jl")
    end

    @require GPUArrays = "0c68f7d7-f131-5f86-a1c3-88cf8149b2d7" begin
      include("../ext/DensityRatioEstimationGPUArraysExt.jl")
    end
  end
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
  KMM,
  uKMM,
  KLIEP,
  LSIF,
  available_optlib,
  default_optlib,
  densratiofunc,
  densratio,

  # estimator fitters
  EstimatorFitter,
  LCV,
  fit

end # module

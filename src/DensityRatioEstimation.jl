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

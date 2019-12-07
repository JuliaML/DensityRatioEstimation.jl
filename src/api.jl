# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DensityRatioEstimator

A method for density ratio estimation.
"""
abstract type DensityRatioEstimator end

"""
    OptimizationLibrary

An optimization library (e.g. Optim.jl, Convex.jl, JuMP.jl).
"""
abstract type OptimizationLibrary end

# optimization libraries for dispatch
struct JuliaLib  <: OptimizationLibrary end
struct OptimLib  <: OptimizationLibrary end
struct ConvexLib <: OptimizationLibrary end
struct JuMPLib   <: OptimizationLibrary end

"""
    density_ratio(x_nu, x_de, dre; [optlib])

Estimate density ratio `p_nu(x) / p_de(x)` with estimator
`dre` and optimization library `optlib`.

Optionally choose an optimization library `optlib` from
the list below:

* `JuliaLib`  - Pure Julia implementation
* `OptimLib`  - Implementation with Optim.jl
* `ConvexLib` - Implementation with Convex.jl
* `JuMPLib`   - Implementation with JuMP.jl
"""
density_ratio(x_nu, x_de, dre::DensityRatioEstimator;
              optlib=_default_optlib(dre)) =
  _density_ratio(x_nu, x_de, dre, optlib)

# internal function with implementation
_density_ratio(x_nu, x_de,
               dre::DensityRatioEstimator,
               opl::Type{O}) where {O<:OptimizationLibrary} =
  @error "not implemented"

# default optimization library for estimator
_default_optlib(dre::DensityRatioEstimator) = @error "not implemented"

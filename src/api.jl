# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DensityRatio

A ratio of two probability density functions.
"""
abstract type DensityRatio end

"""
    evaluate(r, x)

Evaluate density ratio `r` at point `x`.
"""
evaluate(r::DensityRatio, x::AbstractVector) = @error "not implemented"

"""
    DiscreteDensityRatio(points, ratios)

A ratio of two probability density functions for
a discrete set of `points`. `ratios` is a vector
of values for the `points`.
"""
struct DiscreteDensityRatio{P<:AbstractVector,V<:AbstractVector} <: DensityRatio
  points::P
  ratios::V
end

"""
    FunctionalDensityRatio(func)

A ratio of two probability density functions given in
functional form `func`.
"""
struct FunctionalDensityRatio{F<:Function} <: DensityRatio
  func::F
end

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

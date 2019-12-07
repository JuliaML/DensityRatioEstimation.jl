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
    density_ratio(x_nu, x_de, dre; optlib=JuliaLib())

Estimate density ratio `p_nu(x) / p_de(x)` with estimator
`dre` and optimization library `optlib`.
"""
density_ratio(x_nu, x_de, dre::DensityRatioEstimator; optlib=JuliaLib()) =
  _density_ratio(x_nu, x_de, dre, optlib)

# internal function with implementation
_density_ratio(x_nu, x_de,
               dre::DensityRatioEstimator,
               opl::OptimizationLibrary) = @error "not implemented"

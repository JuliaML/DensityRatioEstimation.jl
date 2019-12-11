# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DensityRatioEstimator

A method for density ratio estimation.
"""
abstract type DensityRatioEstimator end

"""
    densratio(x_nu, x_de, dre; [optlib])

Estimate density ratio `p_nu(x) / p_de(x)` with estimator
`dre` and optimization library `optlib`.

Optionally choose an optimization library `optlib` from
the list below:

* `JuliaLib`  - Pure Julia implementation
* `OptimLib`  - Implementation with Optim.jl
* `ConvexLib` - Implementation with Convex.jl
* `JuMPLib`   - Implementation with JuMP.jl
"""
densratio(x_nu, x_de, dre::DensityRatioEstimator;
          optlib=_default_optlib(typeof(dre))) =
  _densratio(x_nu, x_de, dre, optlib)

# internal function with implementation
_densratio(x_nu, x_de, dre::DensityRatioEstimator,
           optlib) = @error "not implemented"

# default optimization library for estimator
_default_optlib(dre::Type{DensityRatioEstimator}) =
  @error "not implemented"

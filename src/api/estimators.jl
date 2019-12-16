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

Estimate density ratio `p_nu(x_de) / p_de(x_de)` with estimator
`dre` and indexable collections of numerator and denominator
samples, `x_nu` and `x_de`.

Optionally choose an optimization library `optlib` from
the list below:

* `JuliaLib`  - Pure Julia implementation
* `OptimLib`  - Implementation with Optim.jl
* `ConvexLib` - Implementation with Convex.jl
* `JuMPLib`   - Implementation with JuMP.jl

See also [`densratiofunc`](@ref).
"""
densratio(x_nu, x_de, dre::DensityRatioEstimator;
          optlib=_default_optlib(typeof(dre))) =
  _densratio(x_nu, x_de, dre, optlib)

"""
    densratiofunc(x_nu, x_de, dre; [optlib])

Similar to [`densratio`](@ref), but return a function
`r(x) = p_nu(x) / p_de(x)` that can be evaluated at a
new unseen sample `x`.

See also [`densratio`](@ref).

### Notes

Only some estimators define a ratio function that can
be evaluated outside `x_de`.
"""
densratiofunc(x_nu, x_de, dre::DensityRatioEstimator;
              optlib=_default_optlib(typeof(dre))) =
  _densratiofunc(x_nu, x_de, dre, optlib)

###################################################
## functions to be implemented by new estimators ##
###################################################

_densratio(x_nu, x_de, dre::DensityRatioEstimator,
           optlib::Type{OptimizationLibrary}) =
  @error "not implemented"

_densratiofunc(x_nu, x_de, dre::DensityRatioEstimator,
               optlib::Type{OptimizationLibrary}) =
  @error "not implemented"

_default_optlib(dre::Type{DensityRatioEstimator}) =
  @error "not implemented"

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
`dre` and optimization library `optlib`. `x_nu` and `x_de`
are indexable collections of numerator and denominator
samples, respectively.

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

"""
    fit(dre, x_nu, x_de; [optlib])

Perform hyperparameter tuning of density ratio
estimator `dre` with numerator and denominator
samples, `x_nu` and `x_de`. Optinally, specify
the optimization library `optlib`.

### Notes

Hyperparameter tuning is not defined for all
density ratio estimators. Therefore, this
function may not work with some estimators.
"""
fit(dre::Type{<:DensityRatioEstimator}, x_nu, x_de;
    optlib=_default_optlib(dre)) =
  _fit(dre, x_nu, x_de, optlib)

###################################################
## functions to be implemented by new estimators ##
###################################################

_densratio(x_nu, x_de, dre::DensityRatioEstimator,
           optlib::Type{OptimizationLibrary}) =
  @error "not implemented"

_fit(dre::Type{DensityRatioEstimator}, x_nu, x_de,
     optlib::Type{OptimizationLibrary}) =
  @error "not implemented"

_default_optlib(dre::Type{DensityRatioEstimator}) =
  @error "not implemented"

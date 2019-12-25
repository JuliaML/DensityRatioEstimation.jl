# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EstimatorFitter

A method to fit a density ratio estimator.
"""
abstract type EstimatorFitter end

"""
    fit(DRE, x_nu, x_de, ranges; [optlib])

Perform hyperparameter tuning of density ratio
estimator `dre` with numerator and denominator
samples, `x_nu` and `x_de` and with hyperparameter
`ranges`. Optionally, specify the optimization
library `optlib`.

### Notes

Hyperparameter tuning is not defined for all
density ratio estimators. Therefore, this
function may not work with some estimators.
"""
fit(dre::Type{<:DensityRatioEstimator}, x_nu, x_de,
    fitter::EstimatorFitter; optlib=_default_optlib(dre)) =
  _fit(dre, x_nu, x_de, fitter, optlib)

################################################
## functions to be implemented by new fitters ##
################################################

_fit(dre::Type{DensityRatioEstimator}, x_nu, x_de,
     fitter::EstimatorFitter, optlib::Type{OptimizationLibrary}) =
  @error "not implemented"

# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DensityRatioEstimator

A method for density ratio estimation.
"""
abstract type DensityRatioEstimator end

density_ratio(x_nu, x_de, dre::DensityRatioEstimator) = @error "not implemented"

# DensityRatioEstimation.jl

[![Build Status](https://travis-ci.com/xukai92/DensityRatioEstimation.jl.svg?branch=master)](https://travis-ci.com/xukai92/DensityRatioEstimation.jl) [![Coverage Status](https://coveralls.io/repos/github/xukai92/DensityRatioEstimation.jl/badge.svg?branch=master)](https://coveralls.io/github/xukai92/DensityRatioEstimation.jl?branch=master)

To get started, read the [introduction](https://htmlpreview.github.io/?https://github.com/xukai92/DensityRatioEstimation.jl/blob/master/docs/intro.html).

## Supported methods

Infinite moment matching based on maximum mean discrepancy (MMD)
- Analytical solution with no constraint
- Numerical solution with optional positivity and normalisation constraints
    - Requires `JuMP.jl` and `Ipopt.jl`.

## References

Sugiyama M, Suzuki T, Kanamori T. Density ratio estimation in machine learning. Cambridge University Press, 2012.

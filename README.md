<p align="center">
  <img src="docs/DensRatio.png" height="200"><br>
  <a href="https://travis-ci.org/JuliaEarth/DensityRatioEstimation.jl">
    <img src="https://travis-ci.org/JuliaEarth/DensityRatioEstimation.jl.svg?branch=master">
  </a>
  <a href="https://codecov.io/gh/JuliaEarth/DensityRatioEstimation.jl">
    <img src="https://codecov.io/gh/JuliaEarth/DensityRatioEstimation.jl/branch/master/graph/badge.svg">
  </a>
  <a href="LICENSE">
    <img src="https://img.shields.io/badge/license-ISC-blue.svg">
  </a>
</p>

Given samples `x_nu` and `x_de` from distributions `p_nu` and `p_de`, it is very
useful to estimate the density ratio `r(x) = p_nu(x) / p_de(x)` for all valid `x`.
This problem is known in the literature as the density ratio estimation problem
([Sugiyama et al. 2012](https://www.cambridge.org/core/books/density-ratio-estimation-in-machine-learning/BCBEA6AEAADD66569B1E85DDDEAA7648)).

Naive solutions based on the ratio of individual estimators for numerator and
denominator densities perform poorly, particularly in high dimensions. This
package provides density ratio estimators that perform well with a moderately
large number of dimensions.

## Installation

Get the latest stable release with Julia's package manager:

```julia
] add DensityRatioEstimation
```

## Usage

Given two indexable collections `x_nu` and `x_de` of samples from `p_nu` and
`p_de`, one can estimate the density ratio at all samples in `x_de`:

```julia
using DensityRatioEstimation, Optim

r = densratio(x_nu, x_de, KLIEP(), optlib=OptimLib)
```

The third argument of the `densratio` function is a density ratio estimator.
Currently, this package implements the following estimators:

| Estimator | References |
| --------- | ---------- |
| KMM       | [Huang et al. 2006](https://papers.nips.cc/paper/3075-correcting-sample-selection-bias-by-unlabeled-data.pdf) |
| KLIEP     | [Sugiyama et al. 2008](https://link.springer.com/article/10.1007/s10463-008-0197-x) |
| LSIF      | [Kanamori et al. 2009](http://www.jmlr.org/papers/volume10/kanamori09a/kanamori09a.pdf) |

The fourth argument `optlib` specifies the optimization package used to implement
the estimator. Some estimators are implemented with different optimization packages
to facilitate the usage in different environments. In the example above, users that
already have the [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl) package in
their environment can promptly use the `KLIEP` estimator implemented with that package.
Each estimator has a default optimization package, and so the function call above
can be simplified given that the optimization package is already loaded:

```julia
r = densratio(x_nu, x_de, KLIEP())
```

Different implementations of the same estimator are loaded using the
[Requires.jl](https://github.com/MikeInnes/Requires.jl) package, and
the keyword argument `optlib` can be any of:

* `JuliaLib`  - Pure Julia implementation
* `OptimLib`  - [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl) implementation
* `ConvexLib` - [Convex.jl](https://github.com/JuliaOpt/Convex.jl) implementation
* `JuMPLib`   - [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl) implementation

To find out the default implementation for an estimator, please use

```julia
default_optlib(KLIEP)
```

and to find out the available implementations, please use

```julia
available_optlib(KLIEP)
```

Some methods support the evaluation of the density ratio at all `x`, besides the
denominator samples. In this case, the following line returns a function `r(x)`
that can be evaluated at new unseen samples:

```julia
r = densratiofunc(x_nu, x_de, KLIEP())
```

### Hyperparameter tuning

Methods like `KLIEP` are equiped with tuning strategies, and its hyperparameters
can be found using the following line:

```julia
dre = fit(KLIEP, x_nu, x_de, LCV((σ=[1.,2.,3.],b=[100]))
```

The function returns a `KLIEP` instance with parameters optimized for the samples.
In this case, the line uses likelihood cross-validation `LCV` as the tuning
strategy. It accepts a named tuple with the hyperparameter ranges for `KLIEP`,
the kernel width `σ` and the number of basis functions `b`. Currently, the
following tuning strategies are implemented:

| Tuning  | References |
| ------- | ---------- |
| LCV     | [Sugiyama et al. 2008](https://link.springer.com/article/10.1007/s10463-008-0197-x) |

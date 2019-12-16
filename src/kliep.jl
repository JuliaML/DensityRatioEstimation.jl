# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    KLIEP(σ=1.0, b=100)

Kullback-Leibler importance estimation procedure (KLIEP).

## Parameters

* `σ` - Bandwidth of Gaussian kernel (default to `1.0`)
* `b` - Number of radial basis functions (default to `100`)

## References

* Sugiyama et al. 2008. Direct importance estimation for
  covariate shift adaptation.

### Author

* Júlio Hoffimann (julio.hoffimann@gmail.com)
"""
@with_kw struct KLIEP{T} <: DensityRatioEstimator
  σ::T=1.0
  b::Int=100
end

_default_optlib(dre::Type{<:KLIEP}) = OptimLib

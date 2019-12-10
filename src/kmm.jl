# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    KMM(σ=1.0, B=1e-4)

Kernel Mean Matching (KMM).

## Parameters

* `σ` - Bandwidth of Gaussian kernel (default to `1.0`)
* `B` - Maximum discrepancy allowed (default to `1e+1`)
* `ϵ` - Tolerance for unit sum (default to `1e-2`)

## References

* Huang et al. 2006. Correcting sample selection bias by
  unlabeled data.
"""
struct KMM{T} <: DensityRatioEstimator
  σ::T
  B::T
  ϵ::T
end

KMM(σ, B) = KMM(σ, B, 1e-2)
KMM(σ) = KMM(σ, 1e+1)
KMM() = KMM(1.0)

_default_optlib(dre::KMM) = JuMPLib

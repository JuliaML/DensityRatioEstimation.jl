# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    KMM(σ=1.0)

Kernel Mean Matching (KMM).

## Parameters

* `σ` - Bandwidth of Gaussian kernel (default to `1.0`)
* `B` - Maximum discrepancy allowed (default to `1.0`)

## References

* Huang et al. 2006. Correcting sample selection bias by
  unlabeled data.
"""
struct KMM{T} <: DensityRatioEstimator
  σ::T
  B::T
end

KMM(σ) = KMM(σ, 1.0)
KMM() = KMM(1.0)

_default_optlib(dre::KMM) = JuMPLib
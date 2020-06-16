# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    KMM(σ=1.0, B=Inf, ϵ=0.01, λ=0.001)

Kernel Mean Matching (KMM).

## Parameters

* `σ` - Bandwidth of Gaussian kernel (default to `2.0`)
* `B` - Maximum discrepancy allowed (default to `Inf`)
* `ϵ` - Tolerance for unit sum (default to `0.01`)
* `λ` - Regularization parameter (default to `0.001`)

## References

* Huang et al. 2006. Correcting sample selection bias by
  unlabeled data.

### Author

* Júlio Hoffimann (julio.hoffimann@gmail.com)
* Kai Xu (xukai921110@gmail.com)
"""
@with_kw struct KMM{T} <: DensityRatioEstimator
  σ::T=2.0
  B::T=Inf
  ϵ::T=0.01
  λ::T=0.001
end

default_optlib(dre::Type{<:KMM}) = JuMPLib

available_optlib(dre::Type{<:KMM}) = [JuliaLib, JuMPLib]

function _kmm_consts(x_nu, x_de, dre::KMM{T}) where {T<:AbstractFloat}
  @unpack σ, λ = dre

  # Gramian matrices for numerator and denominator
  Kdede = gaussian_gramian(x_de; σ=σ)
  if !iszero(λ)
    Kdede += safe_diagm(Kdede, λ)
  end
  Kdenu = gaussian_gramian(x_de, x_nu; σ=σ)

  # number of denominator and numerator samples
  n_de, n_nu = size(Kdenu)

  Kdede, T(n_de / n_nu) * sum(Kdenu, dims=2)
end

function _densratio(x_nu, x_de, dre::KMM,
                    optlib::Type{<:OptimizationLibrary})
  K, κ = _kmm_consts(x_nu, x_de, dre)
  _kmm_ratios(K, κ, dre, optlib)
end

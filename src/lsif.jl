# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LSIF(σ=1.0, b=100, λ=0.0)

Least squares importance fitting.

## Parameters

* `σ` - Bandwidth of Gaussian kernel (default to `1.0`)
* `b` - Maximum number of radial basis functions (default to `100`)
* `λ` - Regularization parameter (default to `0.001`)

## References

* Kanamori et al. 2009. A Least-squares Approach to Direct
Importance Estimation

## Author

* Júlio Hoffimann (julio.hoffimann@gmail.com)
"""
@with_kw struct LSIF{T} <: DensityRatioEstimator
  σ::T=1.0
  b::Int=100
  λ::T=0.001
end

default_optlib(dre::Type{<:LSIF}) = OptimLib

available_optlib(dre::Type{<:LSIF}) = [OptimLib, JuMPLib]

function _densratio(x_nu, x_de, dre::LSIF,
                    optlib::Type{<:OptimizationLibrary})
  c = sample_centers(x_nu, dre.b)
  α = _lsif_coeffs(x_nu, x_de, c, dre, optlib)
  K = gaussian_gramian(x_de, x_nu[c], σ=dre.σ)
  K*α
end

function _densratiofunc(x_nu, x_de, dre::LSIF,
                        optlib::Type{<:OptimizationLibrary})
  c = sample_centers(x_nu, dre.b)
  α = _lsif_coeffs(x_nu, x_de, c, dre, optlib)
  function r(x)
    K = gaussian_gramian([x], x_nu[c], σ=dre.σ)
    dot(K, α)
  end
end

"""
    _lsif_coeffs(x_nu, x_de, basis, dre, optlib)

Return the coefficients of LSIF basis expansion.
"""
_lsif_coeffs(x_nu, x_de, centers::AbstractVector{Int}, dre::LSIF,
             optlib::Type{<:OptimizationLibrary}) =
  @error "not implemented"

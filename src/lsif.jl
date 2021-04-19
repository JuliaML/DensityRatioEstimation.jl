# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LSIF(σ=2.0, b=100, λ=0.001)

Least squares importance fitting.

## Parameters

* `σ` - Bandwidth of Gaussian kernel (default to `2.0`)
* `b` - Maximum number of radial basis functions (default to `10`)
* `λ` - Regularization parameter (default to `0.001`)

## References

* Kanamori et al. 2009. A Least-squares Approach to Direct
Importance Estimation

## Author

* Júlio Hoffimann (julio.hoffimann@gmail.com)
"""
@with_kw struct LSIF{T} <: DensityRatioEstimator
  σ::T=2.0
  b::Int=10
  λ::T=0.001
end

default_optlib(dre::Type{<:LSIF}) = OptimLib

available_optlib(dre::Type{<:LSIF}) = [OptimLib, JuMPLib]

function _densratio(x_nu, x_de, dre::LSIF,
                    optlib::Type{<:OptimizationLibrary})
  K_de, H, h, x_ba = _lsif_consts(x_nu, x_de, dre)
  α = _lsif_coeffs(H, h, dre, optlib)
  K_de*α
end

function _densratiofunc(x_nu, x_de, dre::LSIF,
                        optlib::Type{<:OptimizationLibrary})
  K_de, H, h, x_ba = _lsif_consts(x_nu, x_de, dre)
  α = _lsif_coeffs(H, h, dre, optlib)
  function r(x)
    K = gaussian_gramian([x], x_ba, σ=dre.σ)
    dot(K, α)
  end
end

function _lsif_consts(x_nu, x_de, dre)
  x_ba = sample(x_nu, min(length(x_nu), dre.b), replace=false)
  K_nu = gaussian_gramian(x_nu, x_ba, σ=dre.σ)
  K_de = gaussian_gramian(x_de, x_ba, σ=dre.σ)

  b = length(x_ba)
  Φ = Matrix{eltype(K_de)}(undef, b, b)
  for l′ in 1:b
    φ′ = view(K_de, :, l′)
    for l in 1:l′
      φ = view(K_de, :, l)
      Φ[l,l′] = mean(φ .* φ′)
    end
  end
  H = Symmetric(Φ)
  h = vec(mean(K_nu, dims=1))

  K_de, H, h, x_ba
end

"""
    _lsif_coeffs(x_nu, x_de, basis, dre, optlib)

Return the coefficients of LSIF basis expansion.
"""
_lsif_coeffs(H, h, dre::LSIF, optlib::Type{<:OptimizationLibrary}) =
  @error "not implemented"

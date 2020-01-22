# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    KLIEP(σ=1.0, b=100)

Kullback-Leibler importance estimation procedure (KLIEP).

## Parameters

* `σ` - Bandwidth of Gaussian kernel (default to `1.0`)
* `b` - Maximum number of radial basis functions (default to `100`)

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

default_optlib(dre::Type{<:KLIEP}) = OptimLib

available_optlib(dre::Type{<:KLIEP}) = [OptimLib, ConvexLib]

function _densratio(x_nu, x_de, dre::KLIEP,
                    optlib::Type{<:OptimizationLibrary})
  K_nu, K_de, x_ba = _kliep_consts(x_nu, x_de, dre)
  α = _kliep_coeffs(K_nu, K_de, dre, optlib)
  K_de*α
end

function _densratiofunc(x_nu, x_de, dre::KLIEP,
                        optlib::Type{<:OptimizationLibrary})
  K_nu, K_de, x_ba = _kliep_consts(x_nu, x_de, dre)
  α = _kliep_coeffs(K_nu, K_de, dre, optlib)
  function r(x)
    K = gaussian_gramian([x], x_ba, σ=dre.σ)
    dot(K, α)
  end
end

# constants involved in KLIEP optimization
function _kliep_consts(x_nu, x_de, dre)
  x_ba = sample(x_nu, min(length(x_nu), dre.b), replace=false)
  K_nu = gaussian_gramian(x_nu, x_ba, σ=dre.σ)
  K_de = gaussian_gramian(x_de, x_ba, σ=dre.σ)

  K_nu, K_de, x_ba
end

"""
    _kliep_coeffs(K_nu, K_de, dre, optlib)

Return the coefficients of KLIEP basis expansion.
"""
_kliep_coeffs(K_nu, K_de, dre::KLIEP,
             optlib::Type{<:OptimizationLibrary}) =
  @error "not implemented"

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
  c = _kliep_centers(x_nu, dre)
  α = _kliep_coeffs(x_nu, x_de, c, dre, optlib)
  K = gaussian_gramian(x_de, x_nu[c], σ=dre.σ)
  K*α
end

function _densratiofunc(x_nu, x_de, dre::KLIEP,
                        optlib::Type{<:OptimizationLibrary})
  c = _kliep_centers(x_nu, dre)
  α = _kliep_coeffs(x_nu, x_de, c, dre, optlib)
  function r(x)
    K = gaussian_gramian([x], x_nu[c], σ=dre.σ)
    dot(K, α)
  end
end

"""
    _kliep_centers(x_nu, dre)

Return the indices of `x_nu` used as the kernel centers
in kernel approximation of density ratio function.
"""
function _kliep_centers(x_nu, dre::KLIEP)
  b = dre.b
  n = length(x_nu)
  s = min(n, b)
  sample(1:n, s, replace=false)
end

"""
    _kliep_coeffs(x_nu, x_de, basis, dre, optlib)

Return the coefficients of KLIEP basis expansion.
"""
_kliep_coeffs(x_nu, x_de, centers::AbstractVector{Int}, dre::KLIEP,
             optlib::Type{<:OptimizationLibrary}) =
  @error "not implemented"

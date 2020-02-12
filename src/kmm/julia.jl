# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENSE in the project root.
# ------------------------------------------------------------------

# NOTE: this function is a hack for Zygote.jl compatbility; see lib/zygote.jl
function warn_kmm_julialib(B, ϵ)
  # warn user that closed-form solution does
  # not consider probability simplex constraints
  isinf(B) || @warn "B parameter ignored when optlib=JuliaLib"
  iszero(ϵ) || @warn "ϵ parameter ignored when optlib=JuliaLib"
end

function _densratio(x_nu, x_de, dre::KMM, optlib::Type{JuliaLib})
  # retrieve parameters
  @unpack σ, B, ϵ, λ = dre

  # warn ignored parameters
  warn_kmm_julialib(B, ϵ)

  # Gramian matrices for numerator and denominator
  Kdede = gaussian_gramian(x_de; σ=σ)
  Kdenu = gaussian_gramian(x_de, x_nu; σ=σ)

  # closed-form solution (without constraints)
  _densratio_kmm(Kdede, Kdenu, dre.λ)
end

# Interface to deal with an array of σs for external functions
function _densratio_kmm(Kdede, Kdenu, λ::T) where {T<:AbstractFloat}
  n_de, n_nu = size(Kdenu)
  T(n_de / n_nu) * ((Kdede + safe_diagm(Kdede, λ)) \ vec(sum(Kdenu, dims=2)))
end

# Interface to deal with 
# - CuArrays.jl: slow scalar operations; see lib/cuarrays.jl
# - Zygote.jl: creating a new CuArray is unsupport during pullback; see lib/zygote.jl
safe_diagm(Kdede, λ) = λ * I

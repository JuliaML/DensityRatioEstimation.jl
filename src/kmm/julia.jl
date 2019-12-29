# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function _densratio(x_nu, x_de, dre::KMM, optlib::Type{JuliaLib})
  # retrieve parameters
  @unpack σ, B, ϵ, λ = dre

  # warn user that closed-form solution does
  # not consider probability simplex constraints
  isinf(B) || @warn "B parameter ignored when optlib=JuliaLib"
  iszero(ϵ) || @warn "ϵ parameter ignored when optlib=JuliaLib"

  # number of numerator and denominator samples
  n_nu, n_de = length(x_nu), length(x_de)

  # Gramian matrices for numerator and denominator
  Kdede = gaussian_gramian(x_de, x_de, σ=σ)
  Kdenu = gaussian_gramian(x_de, x_nu, σ=σ)

  # closed-form solution (without constraints)
  (n_de \ n_nu) * (Kdede + λ*I) \ vec(sum(Kdenu, dims=2))
end

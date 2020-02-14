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

function _kmm_ratios(K, κ, dre, optlib::Type{JuliaLib})
  # retrieve parameters
  @unpack B, ϵ = dre

  warn_kmm_julialib(B, ϵ) # warn ignored parameters

  # density ratio via solve
  K \ vec(κ)
end

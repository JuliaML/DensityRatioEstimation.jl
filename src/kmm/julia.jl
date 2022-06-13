# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function _kmm_ratios(K, κ, dre::uKMM, optlib::Type{JuliaLib})
  # density ratio via solver
  K \ vec(κ)
end

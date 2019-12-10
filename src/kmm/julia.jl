# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function _densratio(x_nu, x_de, dre::KMM, optlib::Type{JuliaLib})
  # retrieve parameters
  σ, B, ϵ = dre.σ, dre.B, dre.ϵ

  # number of numerator and denominator samples
  m′, m = length(x_nu), length(x_de)

  # constants for objective and constraints
  K = gaussian_gramian(x_de, x_de, σ=σ)
  A = gaussian_gramian(x_de, x_nu, σ=σ)
  κ = (m / m′) * sum(A, dims=2)

  # TODO: implement closed form solution
end

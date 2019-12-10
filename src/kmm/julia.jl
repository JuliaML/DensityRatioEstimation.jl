# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function _density_ratio(x_nu, x_de, dre::KMM, optlib::Type{JuliaLib})
  # retrieve parameters
  σ = dre.σ

  # number of numerator and denominator samples
  n_nu, n_de = length(x_nu), length(x_de)

  # constants for objective and constraints
  Kdede = gaussian_gramian(x_de, x_de, σ=σ)
  Kdenu = gaussian_gramian(x_de, x_nu, σ=σ)

  # TODO: invert system with \ operator
end

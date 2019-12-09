# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

gaussian_kernel(d, σ) = exp(-d^2/(2σ^2))

abstract type AbstractKMM end

function density_ratio(mmd::AbstractKMM, pdot_dede, pdot_denu, σ)
  Kdede = gaussian_kernel.(pdot_dede, σ)
  Kdenu = gaussian_kernel.(pdot_denu, σ)
  _density_ratio(mmd, Kdede, Kdenu)
end

"""
    density_ratio(mmd::AbstractKMM, x_de, x_nu; σs=nothing)

Estimate `p_nu(x_de) / p_de(x_de)` using mutiple `σ` in `σs`.
If `σs` is not provided, the median of pairwise distances will be used.
"""
function density_ratio(mmd::AbstractKMM, x_de, x_nu; σs=nothing)
  pdot_dede = pairwise(Euclidean(), x_de, dims=2)
  pdot_denu = pairwise(Euclidean(), x_de, x_nu, dims=2)

  if isnothing(σs)
    σ = sqrt(median([pdot_dede..., pdot_denu...]))
    σs = [σ]
  end

  r_de = mapreduce(σ -> density_ratio(mmd, pdot_dede, pdot_denu, σ), +, σs)

  inv(length(σs)) * r_de
end

struct KMMAnalytical{T<:AbstractFloat, S, M} <: AbstractKMM
  ϵ::T
  function KMMAnalytical(ϵ::T=1f-3; method::Symbol=:solve) where {T<:Number}
    @assert method in (:solve, :inv)
    S = iszero(ϵ) ? Val{:false} : Val{:true}
    if T == Int
      ϵ = float(ϵ)
    end
    new{typeof(ϵ), S, Val{method}}(ϵ)
  end
end

function adddiag(mmd::KMMAnalytical{T, Val{:true}, M}, Kdede) where {T, M}
  Kdede + diagm(0 => mmd.ϵ * fill(one(T), size(Kdede, 1)))
end
adddiag(mmd::KMMAnalytical{T, Val{:false}, M}, Kdede) where {T, M} = Kdede

function _density_ratio(mmd::KMMAnalytical{T, S, Val{:inv}}, Kdede, Kdenu) where {T, S}
  n_de, n_nu = size(Kdenu)
  convert(T, n_de / n_nu) * inv(adddiag(mmd, Kdede)) * Kdenu * fill(one(T), n_nu)
end

function _density_ratio(mmd::KMMAnalytical{T, S, Val{:solve}}, Kdede, Kdenu) where {T, S}
  n_de, n_nu = size(Kdenu)
  convert(T, n_de / n_nu) * (adddiag(mmd, Kdede) \ sum(Kdenu; dims=2)[:,1])
end

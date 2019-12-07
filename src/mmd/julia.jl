# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

gaussian_gram_by_pairwise_sqd(pdot, σ) = exp.(-pdot ./ 2(σ ^ 2))
gaussian_gram(x, σ) = gaussian_gram_by_pairwise_sqd(pairwise(Euclidean(), x, dims=2), σ)
gaussian_gram(x, y, σ) = gaussian_gram_by_pairwise_sqd(pairwise(Euclidean(), x, y, dims=2), σ)

abstract type AbstractMMD end

function estimate_ratio(mmd::AbstractMMD, pdot_dede, pdot_denu, σ)
    Kdede = gaussian_gram_by_pairwise_sqd(pdot_dede, σ)
    Kdenu = gaussian_gram_by_pairwise_sqd(pdot_denu, σ)
    return _estimate_ratio(mmd, Kdede, Kdenu)
end

"""
    estimate_ratio(mmd::AbstractMMD, x_de, x_nu; σs=nothing)

Estimate `p_nu(x_de) / p_de(x_de)` using mutiple `σ` in `σs`.
If `σs` is not provided, the median of pairwise distances will be used.
"""
function estimate_ratio(mmd::AbstractMMD, x_de, x_nu; σs=nothing)
    pdot_dede = pairwise(Euclidean(), x_de, dims=2)
    pdot_denu = pairwise(Euclidean(), x_de, x_nu, dims=2)

    if isnothing(σs)
        σ = sqrt(median([pdot_dede..., pdot_denu...]))
        @info "Automatically choose σ using the median of pairwise distances: $σ."
        σs = [σ]
    end

    r_de = mapreduce(σ -> estimate_ratio(mmd, pdot_dede, pdot_denu, σ), +, σs)

    return inv(length(σs)) * r_de
end

struct MMDAnalytical{T<:AbstractFloat, S, M} <: AbstractMMD
    ϵ::T
    function MMDAnalytical(ϵ::T=1f-3; method::Symbol=:solve) where {T<:Number}
        @assert method in (:solve, :inv)
        S = iszero(ϵ) ? Val{:false} : Val{:true}
        if T == Int
            ϵ = float(ϵ)
        end
        return new{typeof(ϵ), S, Val{method}}(ϵ)
    end
end

function adddiag(mmd::MMDAnalytical{T, Val{:true}, M}, Kdede) where {T, M}
    return Kdede + diagm(0 => mmd.ϵ * fill(one(T), size(Kdede, 1)))
end
adddiag(mmd::MMDAnalytical{T, Val{:false}, M}, Kdede) where {T, M} = Kdede

function _estimate_ratio(mmd::MMDAnalytical{T, S, Val{:inv}}, Kdede, Kdenu) where {T, S}
    n_de, n_nu = size(Kdenu)
    return convert(T, n_de / n_nu) * inv(adddiag(mmd, Kdede)) * Kdenu * fill(one(T), n_nu)
end

function _estimate_ratio(mmd::MMDAnalytical{T, S, Val{:solve}}, Kdede, Kdenu) where {T, S}
    n_de, n_nu = size(Kdenu)
    return convert(T, n_de / n_nu) * (adddiag(mmd, Kdede) \ sum(Kdenu; dims=2)[:,1])
end

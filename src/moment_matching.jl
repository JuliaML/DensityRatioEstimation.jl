"""
    pairwise_sqd(x)

Pairwise squared distances between each colum vector in `x`.
"""
function pairwise_sqd(x)
    n = size(x, 2)
    xixj = x' * x
    xsq = sum(x .^ 2; dims=1)
    return repeat(xsq'; outer=(1, n)) + repeat(xsq; outer=(n, 1)) - 2xixj
end

"""
    pairwise_sqd(x, y)

Pairwise squared distances between each colum vector in `x` and `y`.
"""
function pairwise_sqd(x, y)
    nx = size(x, 2)
    ny = size(y, 2)
    xiyj = x' * y
    xsq = sum(x .^ 2; dims=1)
    ysq = sum(y .^ 2; dims=1)
    return repeat(xsq'; outer=(1, ny)) .+ repeat(ysq; outer=(nx, 1)) - 2xiyj
end

gaussian_gram_by_pairwise_sqd(pdot, σ) = exp.(-pdot ./ 2(σ ^ 2))
gaussian_gram(x, σ) = gaussian_gram_by_pairwise_sqd(pairwise_sqd(x), σ)
gaussian_gram(x, y, σ) = gaussian_gram_by_pairwise_sqd(pairwise_sqd(x, y), σ)

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
    pdot_dede = pairwise_sqd(x_de)
    pdot_denu = pairwise_sqd(x_de, x_nu)

    if isnothing(σs)
        σ = sqrt(median([pdot_dede..., pdot_denu...]))
        @info "Automatically choose σ using the median of pairwise distances: $σ."
        σs = [σ]
    end

    r_de = mapreduce(σ -> estimate_ratio(mmd, pdot_dede, pdot_denu, σ), +, σs)
    
    return r_de / convert(Float32, length(σs))
end

@with_kw struct MMDNumerical <: AbstractMMD
    positivity::Bool=true
    normalisation::Bool=true
end

function _estimate_ratio(mmd::MMDNumerical, Kdede, Kdenu)
    n_de, n_nu = size(Kdenu)
    model = Model(with_optimizer(IpoptOptimizer; print_level=0))
    @variable(model, r[1:n_de])
    @objective(model, Min, 1 / n_de ^ 2 * sum(r[i] * Kdede[i,j] * r[j] for i = 1:n_de, j=1:n_de) - 2 / (n_de * n_nu) * sum(r[i] * Kdenu[i,j] for i = 1:n_de, j=1:n_nu))
    if mmd.positivity
        @constraint(model, r .>= 0)
    end
    if mmd.normalisation
        @constraint(model, 1 / n_de * sum(r) == 1)
    end
    optimize!(model)
    return value.(r)
end

struct MMDAnalytical{T<:AbstractFloat,S,M} <: AbstractMMD
    ϵ::T
    function MMDAnalytical(; ϵ::T=1f-3, method::Symbol=:solve) where {T}
        @assert method in (:solve, :inv)
        S = iszero(ϵ) ? Val{:false} : Val{:true}
        if T == Int
            ϵ = float(ϵ)
        end
        new{eltype(ϵ), S, Val{method}}(ϵ)
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

function pairwise_dot(x)
    n = size(x, 2)
    xixj = x' * x
    xsq = sum(x .^ 2; dims=1)
    return repeat(xsq'; outer=(1, n)) + repeat(xsq; outer=(n, 1)) - 2xixj
end

function pairwise_dot(x, y)
    nx = size(x, 2)
    ny = size(y, 2)
    xiyj = x' * y
    xsq = sum(x .^ 2; dims=1)
    ysq = sum(y .^ 2; dims=1)
    return repeat(xsq'; outer=(1, ny)) .+ repeat(ysq; outer=(nx, 1)) - 2xiyj
end

gaussian_gram_by_pairwise_dot(pdot; σ=1) = exp.(-pdot ./ 2(σ ^ 2))
gaussian_gram(x; σ=1) = gaussian_gram_by_pairwise_dot(pairwise_dot(x); σ=σ)
gaussian_gram(x, y; σ=1) = gaussian_gram_by_pairwise_dot(pairwise_dot(x, y); σ=σ)

abstract type AbstractMMD end

function estimate_ratio(mmd::AbstractMMD, pdot_dede, pdot_denu, σ)
    Kdede = gaussian_gram_by_pairwise_dot(pdot_dede; σ=σ)
    Kdenu = gaussian_gram_by_pairwise_dot(pdot_denu; σ=σ)
    return _estimate_ratio(mmd, Kdede, Kdenu)
end

"""
    estimate_ratio(mmd::AbstractMMD, x_de, x_nu; σs=nothing)

Estimate `p_nu(x_de) / p_de(x_de)` using mutiple `σ` in `σs`.
If `σs` is not provided, the median of pairwise distances will be used.
"""
function estimate_ratio(mmd::AbstractMMD, x_de, x_nu; σs=nothing)
    pdot_dede = pairwise_dot(x_de)
    pdot_denu = pairwise_dot(x_de, x_nu)

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

@with_kw struct MMDAnalytical{T} <: AbstractMMD
    ϵ::T=1f-3
end

function _estimate_ratio(mmd::MMDAnalytical{T}, Kdede, Kdenu) where {T}
    n_de, n_nu = size(Kdenu)
    return convert(T, n_de / n_nu) * inv(Kdede + diagm(0 => mmd.ϵ * ones(T, n_de))) * Kdenu * ones(T, n_nu) 
end

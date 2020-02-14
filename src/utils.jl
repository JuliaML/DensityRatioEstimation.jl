# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    euclidsq(x, y)

Return the Euclidean distance between two indexable objects.
"""
euclidsq(x, y) = sum((x[i] - y[i])^2 for i in eachindex(x))

# Support matrix data in a GPU and AD compatible way

function euclidsq(X::T) where {T<:AbstractMatrix}
    n = size(X, 2)
    xixj = transpose(X) * X
    xsq = sum(X .^ 2; dims=1)
    transpose(xsq) .+ xsq - 2xixj
end

function euclidsq(X::T, Y::T) where {T<:AbstractMatrix}
    nx = size(X, 2)
    ny = size(Y, 2)
    xiyj = transpose(X) * Y
    xsq = sum(X .^ 2; dims=1)
    ysq = sum(Y .^ 2; dims=1)
    transpose(xsq) .+ ysq - 2xiyj
end

"""
    gaussian_gramian(xs, ys, σ=1)

Gramian matrix for samples `xs` and `ys` using a Gaussian kernel
kernel with bandwidth `σ`.
"""
gaussian_gramian(xs; kwargs...) = gaussian_gramian(xs, xs; kwargs...)
gaussian_gramian(xs, ys; σ=1) =
  [exp(-euclidsq(x, y) / 2σ^2) for x in xs, y in ys]

function gaussian_gramian(X::T, Y::T; σ=1) where {T<:AbstractMatrix}
    gaussian_gramian_by_euclidsq(euclidsq(X, Y), σ)
end

gaussian_gramian_by_euclidsq(esq, σ) = exp.(-esq ./ 2σ^2)

"""
    safe_diagm(mat, a)

Generate a squared matrix whose diagonal is `a` that is 
compatible to perform addition on `mat`. It hebaves 
differently based on `mat` is on CPU or GPU.

It is compatible with
- CuArrays.jl (see lib/cuarrays.jl)
- Zygote.jl (see lib/zygote.jl)
"""
safe_diagm(mat, a) = a * I

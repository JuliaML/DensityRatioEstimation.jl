# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    euclidsq(x, y)

Return the Euclidean distance between two indexable objects.
"""
euclidsq(x, y) = sum((x[i] - y[i])^2 for i in eachindex(x))

# Support matrix data in a GPU and AD compatible way

"""
    euclidsq(X::T, Y::T) where {T<:AbstractMatrix}

Return the Euclidean distance between columns in two matrices.
"""
function euclidsq(X::T, Y::T) where {T<:AbstractMatrix}
  XiXj = transpose(X) * Y
  x² = sum(X .^ 2; dims=1)
  y² = sum(Y .^ 2; dims=1)
  transpose(x²) .+ y² - 2XiXj
end

"""
    euclidsq(X::T) where {T<:AbstractMatrix}

Effective version of `euclidsq(X, X)`.
"""
function euclidsq(X::T) where {T<:AbstractMatrix}
  XiXj = transpose(X) * X
  x² = sum(X .^ 2; dims=1)
  transpose(x²) .+ x² - 2XiXj
end

"""
    gaussian_gramian(xs, ys, σ=1)

Gramian matrix for samples `xs` and `ys` using a Gaussian kernel
kernel with bandwidth `σ`.
"""
gaussian_gramian(xs; kwargs...) = gaussian_gramian(xs, xs; kwargs...)
gaussian_gramian(xs, ys; σ=1) = [exp(-euclidsq(x, y) / 2σ^2) for x in xs, y in ys]

function gaussian_gramian(X::T, Y::T; σ=1) where {T<:AbstractMatrix}
  gaussian_gramian(euclidsq(X, Y), σ)
end

gaussian_gramian(esq, σ::AbstractFloat) = exp.(-esq ./ 2σ^2)

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

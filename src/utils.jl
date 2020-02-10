# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    euclidsq(x, y)

Return the Euclidean distance between two indexable objects.
"""
euclidsq(x, y) = sum((x[i] - y[i])^2 for i in eachindex(x))

"""
    gaussian_gramian(xs, ys, σ=1)

Gramian matrix for samples `xs` and `ys` using a Gaussian kernel
kernel with bandwidth `σ`.
"""
gaussian_gramian(xs; kwargs...) = gaussian_gramian(xs, xs; kwargs...)
gaussian_gramian(xs, ys; σ=1) =
  [exp(-euclidsq(x, y) / 2σ^2) for x in xs, y in ys]
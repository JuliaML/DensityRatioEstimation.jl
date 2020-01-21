# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
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
function gaussian_gramian(xs, ys; σ=1)
  [exp(-euclidsq(x, y) / 2σ^2) for x in xs, y in ys]
end

"""
    select_centers(x_nu, nmax)

Return the indices of `x_nu` used as the kernel centers
in kernel approximation of density ratio function. The
parameter `nmax` is the maximum number of indices to select.
"""
function select_centers(x_nu, nmax::Int)
  n = length(x_nu)
  s = min(n, nmax)
  sample(1:n, s, replace=false)
end

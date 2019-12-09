# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function gaussian_gramian(xs, ys; σ=1)
  X = reduce(hcat, xs); Y = reduce(hcat, ys)
  D = pairwise(Euclidean(), X, Y, dims=2)
  @. exp(-D^2/2σ^2)
end

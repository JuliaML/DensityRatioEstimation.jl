# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

euclidsq(x, y) = sum((x[i] - y[i])^2 for i in eachindex(x))

function gaussian_gramian(xs, ys; σ=1)
  [exp(-euclidsq(x, y) / 2σ^2) for x in xs, y in ys]
end

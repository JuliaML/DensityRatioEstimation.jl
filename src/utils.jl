# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function gaussian_gramian(xs, ys; σ=1)
	euclidsq(x, y) = sum((x .- y).^2)
	return [exp(-euclidsq(x, y) / 2σ^2) for x in xs, y in ys]
end

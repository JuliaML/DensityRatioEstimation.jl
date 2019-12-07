# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DensityRatio

A ratio of two probability density functions.
"""
abstract type DensityRatio end

"""
    evaluate(r, x)

Evaluate density ratio `r` at point `x`.
"""
evaluate(r::DensityRatio, x::AbstractVector) = @error "not implemented"

"""
    DiscreteDensityRatio(points, ratios)

A ratio of two probability density functions for
a discrete set of `points`. `ratios` is a vector
of values for the `points`.
"""
struct DiscreteDensityRatio{P<:AbstractVector,V<:AbstractVector} <: DensityRatio
  points::P
  ratios::V
end

"""
    FunctionalDensityRatio(func)

A ratio of two probability density functions given in
functional form `func`.
"""
struct FunctionalDensityRatio{F<:Function} <: DensityRatio
  func::F
end

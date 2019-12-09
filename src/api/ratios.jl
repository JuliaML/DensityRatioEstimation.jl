# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DensityRatio

A ratio of two probability density functions.
"""
abstract type DensityRatio end

#------------------
# IMPLEMENTATIONS
#------------------
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

# TODO: implement interpolation scheme for discrete density ratio
(r::DiscreteDensityRatio)(x) = @error "not implemented"

"""
    getindex(r, inds)

Return the ratio values stored in the indices
`inds` of discrete density ratio `r`.
"""
Base.getindex(r::DiscreteDensityRatio, inds) = r.ratios[inds]

"""
    FunctionalDensityRatio(func)

A ratio of two probability density functions given in
functional form `func`.
"""
struct FunctionalDensityRatio{F<:Function} <: DensityRatio
  func::F
end

(r::FunctionalDensityRatio)(x) = r.func(x)

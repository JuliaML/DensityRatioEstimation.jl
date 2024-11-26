# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    euclidsq(x, y)

Return the Euclidean distance between two indexable objects.
"""
euclidsq(x, y) = sum((x[i] - y[i])^2 for i in eachindex(x))

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
compatible to perform addition on `mat`. It behaves 
differently based on whether `mat` is on a CPU or GPU.
"""
safe_diagm(mat, a) = a * I

# avoid `mat + a * I` on GPU which involves scalar operations and is slow
function safe_diagm(mat::AbstractGPUMatrix, a)
  diag = similar(mat, size(m, 1))
  fill!(diag, a)
  Diagonal(diag)
end

@non_differentiable safe_diagm(::Any, ::Any)

###################################################
##   Functions and objects for throwing errors   ##
###################################################

OPTLIB_DICT = Dict("JuliaLib" => "Julia", "OptimLib" => "Optim", "ConvexLib" => "Convex", "JuMPLib" => "JuMP")

function _throw_opt_error(dre::DensityRatioEstimator, optlib::Type{<:OptimizationLibrary})
  dre_name = nameof(typeof(dre))
  lib_name = OPTLIB_DICT[string(optlib)]
  optlib_options = join(available_optlib(dre), ", ")
  error(
    "Attempted to compute $(dre_name) density ratios using (possibly default) optimization library $(optlib), but this library has either not been loaded or is not implemented for use with $(dre_name). Available options for `optlib`: $(optlib_options). If $(optlib) is contained within the available options, be sure to call `using $(lib_name)` before calling `densratio`, `densratiofunc`, or other functions fitting the $(dre_name) estimator."
  )
end

function _throw_not_implemented_error(func::String, dre::Type{<:DensityRatioEstimator})
  dre_name = nameof(dre)
  error(
    "Attempted to call `$(func)($(dre_name))` but this function has not been implemented for density ratio estimator of type $(dre_name)."
  )
end

_throw_not_fit_error(dre::Type{<:DensityRatioEstimator}, fitter::EstimatorFitter, optlib::Type{<:OptimizationLibrary}) =
  error(
    "Attempted to `fit` estimator `$(nameof(dre))` using fitter `$(typeof(fitter))` with optimization library `optlib=$(optlib)`, but no `fit` function has been implemented for this combination of estimator, fitter, and optimization library."
  )

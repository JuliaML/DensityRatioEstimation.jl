# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module DensityRatioEstimationGPUArraysExt

using DensityRatioEstimation
using GPUArrays
using LinearAlgebra

# Avoid `mat + a * I` with CUDA which involves scalar operations and is slow
function DensityRatioEstimation.safe_diagm(mat::M, a::T) where {M<:GPUArrays.AbstractGPUArray{T,2}} where {T}
  diag = similar(mat, size(m, 1))
  fill!(diag, a)
  LinearAlgebra.Diagonal(diag)
end

end #module

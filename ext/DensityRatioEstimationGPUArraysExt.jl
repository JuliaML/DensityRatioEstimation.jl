# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------
module DensityRatioEstimationGPUArraysExt

    if isdefined(Base, :get_extension)
        using DensityRatioEstimation
        using GPUArrays
    else
        using ..DensityRatioEstimation
        using ..GPUArrays
    end
    using LinearAlgebra

    # Aviod `mat + a * I` with CUDA which involes scalar operations and is slow
    function DensityRatioEstimation.safe_diagm(mat::<:M, a::T) where {M<:GPUArrays.AbstractGPUArray{T, 2},T}
        LinearAlgebra.Diagonal(M(fill(a,size(mat,1))))
    end

end #module

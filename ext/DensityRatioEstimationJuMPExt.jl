# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------
module DensityRatioEstimationJuMPExt

    if isdefined(Base, :get_extension)
        using DensityRatioEstimation
        using DensityRatioEstimation: LSIF, JuMPLib, AbstractKMM, uKMM, KMM
        using DensityRatioEstimation.Parameters
        using JuMP
        using Ipopt
    else
        using ..DensityRatioEstimation
        using ..DensityRatioEstimation: LSIF, JuMPLib, AbstractKMM, uKMM, KMM
        using ..DensityRatioEstimation.Parameters
        using ..JuMP
        using ..Ipopt
    end
    using LinearAlgebra
    
    include("../src/kmm/jump.jl")
    include("../src/lsif/jump.jl")

end #module

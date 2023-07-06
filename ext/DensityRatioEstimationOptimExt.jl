# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module DensityRatioEstimationOptimExt

using DensityRatioEstimation
using DensityRatioEstimation: KLIEP, LSIF, OptimLib
using Optim

using LinearAlgebra

include("../src/kliep/optim.jl")
include("../src/lsif/optim.jl")

end #module

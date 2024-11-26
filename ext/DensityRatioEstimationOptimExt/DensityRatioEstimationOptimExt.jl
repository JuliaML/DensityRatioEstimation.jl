# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module DensityRatioEstimationOptimExt

using DensityRatioEstimation: KLIEP, LSIF, OptimLib
using LinearAlgebra
using Optim

import DensityRatioEstimation

include("kliep.jl")
include("lsif.jl")

end #module

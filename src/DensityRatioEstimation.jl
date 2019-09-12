module DensityRatioEstimation

using Statistics, LinearAlgebra, JuMP
import Ipopt; IpoptOptimizer = Ipopt.Optimizer

include("moment_matching.jl")
export estimate_r_mmd, get_r_hat_numerically, get_r_hat_analytical

end

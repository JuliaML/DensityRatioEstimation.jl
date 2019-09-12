module DensityRatioEstimation

using Parameters, Statistics, LinearAlgebra, JuMP
import Ipopt; IpoptOptimizer = Ipopt.Optimizer

include("moment_matching.jl")
export MMDAnalytical, MMDNumerical

export estimate_ratio

end

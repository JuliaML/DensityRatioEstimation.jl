module DensityRatioEstimation

using Parameters, Statistics, LinearAlgebra, JuMP
import Ipopt; IpoptOptimizer = Ipopt.Optimizer

estimate_ratio() = error("`estimate_ratio()` is not implemented.")
export estimate_ratio

include("moment_matching.jl")
export MMDAnalytical, MMDNumerical

end

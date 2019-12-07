using DensityRatioEstimation
using Distributions
using Optim
using JuMP, Ipopt
using Test

# environment settings
islinux = Sys.islinux()
istravis = "TRAVIS" ∈ keys(ENV)
datadir = joinpath(@__DIR__,"data")

# list of tests
testfiles = [
  "mmd.jl",
  "kliep.jl"
]

@testset "DensityRatioEstimation.jl" begin
  for testfile in testfiles
    include(testfile)
  end
end

using DensityRatioEstimation
using Distributions
using LinearAlgebra
using Statistics
using Optim
using JuMP, Ipopt
using Convex, ECOS
using Plots, VisualRegressionTests
using Test, Pkg, Random

# workaround GR warnings
ENV["GKSwstype"] = "100"

# environment settings
islinux = Sys.islinux()
istravis = "TRAVIS" ∈ keys(ENV)
datadir = joinpath(@__DIR__,"data")
visualtests = !istravis || (istravis && islinux)
if !istravis
  Pkg.add("Gtk")
  using Gtk
end

# helper funcions
include("utils.jl")

# simple cases for testing
pair₁ = Normal(1,1), Normal(0,2)
pair₂ = MixtureModel([Normal(-2,1),Normal(2,2)], [0.2,0.8]), Normal(0,2)

# list of tests
testfiles = [
  "basic.jl",
  "kmm.jl",
  "kliep.jl",
  "lsif.jl"
]

@testset "DensityRatioEstimation.jl" begin
  for testfile in testfiles
    include(testfile)
  end
end

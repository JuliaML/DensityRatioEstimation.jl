using DensityRatioEstimation
using Distributions
using LinearAlgebra
using Statistics
using Test, StableRNGs
using ReferenceTests, ImageIO
import CairoMakie as Mke

@testset "Throw errors: optlib not loaded" begin

  @test_throws ErrorException densratio([0], [0], KLIEP())
  @test_throws ErrorException densratiofunc([0], [0], KLIEP())
  @test_throws ErrorException densratio([0], [0], LSIF())
  @test_throws ErrorException densratiofunc([0], [0], LSIF())
  @test_throws ErrorException densratio([0], [0], KMM(); optlib = JuMPLib)

end

using Optim
using JuMP, Ipopt
using Convex, ECOS

@testset "Throw errors: optlib not available" begin

  @test_throws ErrorException densratio([0], [0], KLIEP(); optlib = JuMPLib)
  @test_throws ErrorException densratio([0], [0], LSIF(); optlib = ConvexLib)
  @test_throws ErrorException densratio([0], [0], KMM(); optlib = ConvexLib)

end

@testset "Throw errors: undefined functions" begin

  @test_throws ErrorException fit(KMM, [0], [0], LCV((σ=[1.],b=[10])))

  struct NewDRE <: DensityRatioEstimator end

  @test_throws ErrorException densratio([0], [0], NewDRE())
  @test_throws ErrorException densratiofunc([0], [0], NewDRE())
  @test_throws ErrorException default_optlib(NewDRE)
  @test_throws ErrorException available_optlib(NewDRE)

end

# environment settings
isCI = "CI" ∈ keys(ENV)
islinux = Sys.islinux()
visualtests = !isCI || (isCI && islinux)
datadir = joinpath(@__DIR__, "data")

# helper functions
include("utils.jl")

# simple cases for testing
pair₁ = Normal(1, 1), Normal(0, 2)
pair₂ = MixtureModel([Normal(-2, 1), Normal(2, 2)], [0.2, 0.8]), Normal(0, 2)

# list of tests
testfiles = ["basic.jl", "kmm.jl", "kliep.jl", "lsif.jl"]

@testset "DensityRatioEstimation.jl" begin
  for testfile in testfiles
    include(testfile)
  end
end

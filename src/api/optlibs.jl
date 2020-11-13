# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    OptimizationLibrary

An optimization library (e.g. Optim.jl, Convex.jl, JuMP.jl).
"""
abstract type OptimizationLibrary end

# optimization libraries for dispatch
struct JuliaLib  <: OptimizationLibrary end
struct OptimLib  <: OptimizationLibrary end
struct ConvexLib <: OptimizationLibrary end
struct JuMPLib   <: OptimizationLibrary end

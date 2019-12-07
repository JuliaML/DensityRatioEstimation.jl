using Test

@testset "Tests" begin
    tests = [
        "mmd.jl",
    ]

    for t in tests
        include(t)
    end
end

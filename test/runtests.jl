using Test

@testset "Tests" begin
    tests = [
        "moment_matching.jl",
    ]

    for t in tests
        include(t)
    end
end

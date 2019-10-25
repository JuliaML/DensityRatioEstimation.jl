module DensityRatioEstimation

using Requires, Parameters, Statistics, LinearAlgebra

include("moment_matching.jl")
export MMDAnalytical, MMDNumerical

export estimate_ratio

function __init__()
    # The path name `glue` comes from here: https://github.com/JuliaLang/Pkg.jl/issues/1285#issuecomment-525891481
    @require JuMP="4076af6c-e467-56ae-b986-b466b2749572" @require Ipopt="b6b21f68-93f8-5de0-b562-5493be1d77c9" include("glue/opt.jl")
end

end

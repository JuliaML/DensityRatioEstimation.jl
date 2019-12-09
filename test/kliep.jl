# @testset "KLIEP" begin
#   @testset "Optim" begin
#
#   end
# end

using Revise
using DensityRatioEstimation
using Optim
using Convex, ECOS
using JuMP, Ipopt

using Distributions
using Plots
using Random

Random.seed!(123)

n_nu = 1500
n_de = 1000
# d_nu = Normal(0,1)
d_nu = MixtureModel([Normal(-2,1), Normal(2,2)], [0.2, 0.8])
d_de = Normal(0,2)
x_nu = rand(d_nu, n_nu)
x_de = rand(d_de, n_de)

# analytical density ratio
r = pdf.(d_nu, x_de) ./ pdf.(d_de, x_de)

# density ratio estimate
# r̂ = @time density_ratio(x_nu, x_de, KLIEP(), optlib=ConvexLib)
r̂ = @time density_ratio(x_nu, x_de, KMM(), optlib=JuMPLib)

# sort samples for plotting
ind = sortperm(x_de)
xs = x_de[ind]
rs = r[ind]
r̂s = r̂[ind]
# r̂s = [r̂(x) for x in xs]

ps = pdf.(d_nu, xs)
p̂s = r̂s .* pdf.(d_de, xs)

plot(xs, ps)
plot!(xs, p̂s)

# @testset "KLIEP" begin
#   @testset "Optim" begin
#
#   end
# end
# using Plots

n_nu = 1500
n_de = 1000
d_nu = Normal(0,1)
d_de = Normal(0,2)
x_nu = rand(d_nu, n_nu)
x_de = rand(d_de, n_de)

# analytical density ratio
r = pdf.(d_nu, x_de) ./ pdf.(d_de, x_de)

# density ratio estimate
r̂ = density_ratio(x_nu, x_de, KLIEP())

# sort samples for plotting
ind = sortperm(x_de)
xs = x_de[ind]
rs = r[ind]
r̂s = r̂[ind]

ps = pdf.(d_nu, xs)
p̂s = r̂s.* pdf.(d_de, xs)

# plot(xs, ps)
# plot!(xs, p̂s)

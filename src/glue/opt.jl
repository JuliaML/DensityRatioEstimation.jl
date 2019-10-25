using JuMP
import Ipopt; IpoptOptimizer = Ipopt.Optimizer

@with_kw struct MMDNumerical <: AbstractMMD
    positivity::Bool=true
    normalisation::Bool=true
end

function _estimate_ratio(mmd::MMDNumerical, Kdede, Kdenu)
    n_de, n_nu = size(Kdenu)
    model = Model(with_optimizer(IpoptOptimizer; print_level=0))
    @variable(model, r[1:n_de])
    @objective(model, Min, 1 / n_de ^ 2 * sum(r[i] * Kdede[i,j] * r[j] for i = 1:n_de, j=1:n_de) - 2 / (n_de * n_nu) * sum(r[i] * Kdenu[i,j] for i = 1:n_de, j=1:n_nu))
    if mmd.positivity
        @constraint(model, r .>= 0)
    end
    if mmd.normalisation
        @constraint(model, 1 / n_de * sum(r) == 1)
    end
    optimize!(model)
    return value.(r)
end
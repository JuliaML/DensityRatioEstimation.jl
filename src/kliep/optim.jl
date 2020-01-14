# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

using .Optim

function _kliep_coeffs(x_nu, x_de, centers::AbstractVector{Int},
                       dre::KLIEP, optlib::Type{OptimLib})
  # retrieve parameters
  σ, b = dre.σ, length(centers)

  # number of numerator and denominator samples
  n_nu, n_de = length(x_nu), length(x_de)

  # constants for objective
  K = gaussian_gramian(x_nu, x_nu[centers], σ=σ)

  # constants for equality constraints
  P = gaussian_gramian(x_de, x_nu[centers], σ=σ)
  A = sum(P, dims=1)
  lc = uc = [n_de]

  # constants for inequality constraints
  T = eltype(x_de[1])
  lx = fill(zero(T), b)
  ux = fill(Inf, b)

  # objective
  f(α) = -sum(log, K*α)
  function ∇f!(g, α)
    p = K*α
    for l in 1:b
      g[l] = -sum(K[j,l] / p[j] for j in 1:n_nu)
    end
  end
  function ∇²f!(h, α)
    p = K*α
    for k in 1:b, l in 1:b
      h[k,l] = sum(view(K,:,k) .* view(K,:,l) ./ p)
    end
  end

  # equality constraint
  c!(c, α)    = c  .= A*α
  J!(J, α)    = J  .= A
  H!(H, α, λ) = H .+= 0

  # initial guess
  αₒ = fill(n_de/sum(A), b)

  # optimization problem
  objective   = TwiceDifferentiable(f, ∇f!, ∇²f!, αₒ)
  constraints = TwiceDifferentiableConstraints(c!, J!, H!, lx, ux, lc, uc)
  initguess   = αₒ

  # solve problem with interior-point primal-dual Newton
  solution = optimize(objective, constraints, initguess, IPNewton())

  # optimal coefficients
  solution.minimizer
end

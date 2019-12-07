# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

using .Optim

function _density_ratio(x_nu, x_de, dre::KLIEP, opl::OptimLib)
  # retrieve parameters
  σ, b = dre.σ, dre.b

  # Gaussian kernel
  kern(x, y) = exp(-norm(x - y)^2 / (2σ^2))

  # number of numerator and denominator samples
  n_nu = length(x_nu)
  n_de = length(x_de)

  # basis for kernel approximation
  basis = sample(1:n_nu, b, replace=false)

  # constants for objective
  Φ = Matrix{Float64}(undef, n_nu, b)
  for l in 1:b
    xₗ = x_nu[basis[l]]
    for k in 1:n_nu
      xₖ = x_nu[k]
      Φ[k,l] = kern(xₖ, xₗ)
    end
  end

  # constants for equality constraints
  A = Matrix{Float64}(undef, 1, b)
  for l in 1:b
    xₗ = x_nu[basis[l]]
    A[l] = sum(kern(x_de[k], xₗ) for k in 1:n_de)
  end
  lc = uc = [n_de]

  # constants for inequality constraints
  lx = fill(0.0, b)
  ux = fill(Inf, b)

  # objective
  f(α) = -sum(log, Φ*α)
  function ∇f!(g, α)
    p = Φ*α
    for l in 1:b
      g[l] = -sum(Φ[j,l] / p[j] for j in 1:n_nu)
    end
  end
  function ∇²f!(h, α)
    p = Φ*α
    for k in 1:b, l in 1:b
      h[k,l] = sum(view(Φ,:,k) .* view(Φ,:,l) ./ p)
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

  # optimal weights
  α = solution.minimizer
  weights = map(1:n_de) do i
    sum(α[l] * kern(x_de[i], x_nu[c]) for (l, c) in enumerate(basis))
  end

  weights
end

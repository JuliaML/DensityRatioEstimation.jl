# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENSE in the project root.
# ------------------------------------------------------------------

using .Optim

function _lsif_coeffs(H, h, dre::LSIF, optlib::Type{OptimLib})
  # retrieve parameters
  λ, b = dre.λ, length(h)

  # constants for inequality constraints
  T  = eltype(H)
  lx = fill(zero(T), b)
  ux = fill(Inf, b)

  # objective
  f(α)       = (1/2) * dot(α, H*α - 2h) + λ * sum(α)
  ∇f!(grad, α)  = (grad .= H*α .- h .+ λ)
  ∇²f!(hess, α) = (hess .= H)

  # initial guess
  αₒ = fill(one(T), b)

  # optimization problem
  objective   = TwiceDifferentiable(f, ∇f!, ∇²f!, αₒ)
  constraints = TwiceDifferentiableConstraints(lx, ux)
  initguess   = αₒ

  # solve problem with interior-point primal-dual Newton
  solution = optimize(objective, constraints, initguess, IPNewton())

  # optimal coefficients
  solution.minimizer
end

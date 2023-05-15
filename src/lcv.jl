# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LCV(ranges, nfolds=10)

Likelihood cross-validation with parameter ranges
`ranges` and number of folds `nfolds`.

## References

* Sugiyama et al. 2008. Direct importance estimation for
  covariate shift adaptation.

### Author

* Júlio Hoffimann (julio.hoffimann@gmail.com)
"""
struct LCV <: EstimatorFitter
  ranges::NamedTuple
  nfolds::Int
end

LCV(ranges::NamedTuple) = LCV(ranges, 10)

function _fit(::Type{<:KLIEP}, x_nu, x_de, fitter::EstimatorFitter, optlib::Type{<:OptimizationLibrary})
  # retrieve parameters
  ranges = fitter.ranges
  nfolds = fitter.nfolds
  npts = length(x_nu)

  @assert nfolds ≤ npts "number of folds must be smaller than number of numerator samples"

  # partition numerator samples into folds
  folds = collect(Iterators.partition(1:npts, npts ÷ nfolds))

  # initialize loss and optimal estimator
  Ĵmax, dre⭐ = -Inf, KLIEP()

  # loop over hyperparameters
  for σ in ranges.σ, b in ranges.b
    # density ratio estimator
    dre = KLIEP(σ=σ, b=b)

    # estimate loss with cross-validation
    Ĵₖ = map(1:nfolds) do k
      # training and hold-out samples
      train = [ind for i in vcat(1:(k - 1), (k + 1):nfolds) for ind in folds[i]]
      hold = folds[k]

      # perform estimation with training samples
      r = densratiofunc(x_nu[train], x_de, dre, optlib=optlib)

      # evaluate loss with hold-out samples
      mean(log(r(x_nu[j])) for j in hold)
    end

    # mean loss for hyperparameters
    Ĵ = mean(Ĵₖ)

    # update and continue
    if Ĵ > Ĵmax
      Ĵmax = Ĵ
      dre⭐ = dre
    end
  end

  # optimal estimator
  dre⭐
end

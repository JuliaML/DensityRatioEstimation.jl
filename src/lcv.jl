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

function _fit(::Type{<:KLIEP}, x_nu, x_de,
              fitter::EstimatorFitter,
              optlib::Type{<:OptimizationLibrary})
  # retrieve parameters
  ranges = fitter.ranges
  nfolds = fitter.nfolds
  n_nu   = length(x_nu)

  @assert nfolds ≤ n_nu "number of folds must be smaller than number of numerator samples"

  # partition numerator samples into folds
  folds = collect(Iterators.partition(1:n_nu, n_nu ÷ nfolds))

  # loop over hyperparameters
  for σ in ranges.σ, b in ranges.b
    # estimate loss with cross-validation
    for k in 1:nfolds
      # training and hold-out samples
      train = [ind for i in vcat(1:k-1, k+1:nfolds) for ind in folds[i]]
      hold  = folds[k]

      # perform density ratio estimation with training samples
      r = densratiofunc(x_nu[train], x_de, KLIEP(σ=σ, b=b), optlib=optlib)

      # evaluate ratio function at hold-out samples
      w = [r(x_nu[j]) for j in hold]
    end
  end
end

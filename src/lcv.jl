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
  ranges = fitter.ranges
  nfolds = fitter.nfolds
  npts   = length(x_nu)

  @assert nfolds ≤ npts "number of folds must be smaller than number of points"

  folds = collect(Iterators.partition(1:npts, npts ÷ nfolds))

  for σ in ranges.σ, b in ranges.b
    for k in 1:nfolds
      hold  = folds[k]
      train = folds[vcat(1:k-1,k+1:nfolds)]
      r = _densratio(x_nu[train], x_de, KLIEP(σ=σ, b=b), optlib)
    end
  end
end

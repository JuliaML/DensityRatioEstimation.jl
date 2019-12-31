function plot_d_nu(pair, x_de, r̂)
  # extract pair
  d_nu, d_de = pair

  # sort vectors for plotting
  inds = sortperm(x_de)
  xs = x_de[inds]
  r̂s = r̂[inds]

  # actual and estimated d_nu
  ps = pdf.(d_nu, xs)
  p̂s = r̂s .* pdf.(d_de, xs)

  plot(xs, ps, label="d_nu (actual)")
  plot!(xs, p̂s, label="d_nu (estimate)")
end

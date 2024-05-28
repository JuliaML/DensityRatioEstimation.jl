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

  fig = Mke.Figure(size=(800, 800))
  ax = Mke.Axis(fig[1, 1])
  Mke.lines!(ax, xs, ps, label="d_nu (actual)")
  Mke.lines!(ax, xs, p̂s, label="d_nu (estimate)")
  Mke.axislegend(ax)
  fig
end

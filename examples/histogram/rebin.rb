#!/usr/bin/env ruby
require("gsl")

h = GSL::Histogram.alloc(100, [-4, 4])

rng = GSL::Rng.alloc
for i in 0..10000 do
  h.increment(rng.gaussian)
end

h2 = h.rebin
h5 = h.rebin(5)
h7 = h.rebin(7)

printf("%d %d %d %d\n", h.n, h2.n, h5.n, h7.n)
printf("%f %f %f %f\n", h.sigma, h2.sigma, h5.sigma, h7.sigma)
GSL::graph(h, h2, h5, h7, "-T X -C -g 3 -x -4 4")

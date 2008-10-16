#!/usr/bin/env ruby
require("gsl")

N = 10000
BINS = 100

rng1 = GSL::Rng.alloc(2)

h = GSL::Histogram.alloc(BINS, [-5, 5])

for i in 0...N do
  r1 = rng1.gaussian
  h.increment(r1)
end

# Integrate: cumulative distribution
hi = h.integrate

a = hi.diff

# Scale the histograms to ~ 1 at the maximum (to display together)
h.scale!(1.0/h[BINS/2])
a.scale!(1.0/a[BINS/2])

hi.normalize!                  # this is equivalent to hi.scale(1.0/hi[BINS-1])

GSL::graph(h, hi,a,  "-T X -C -g 3")


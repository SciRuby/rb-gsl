#!/usr/bin/env ruby
require("gsl")
rng = GSL::Rng.alloc()

N = 100
XMAX = 5.0
sigma = 1.0
h = GSL::Histogram.alloc(N, 0, XMAX)
for i in 0...10000 do
  x = rng.rayleigh(sigma)
  h.increment(x)
end

sig, amp = h.fit_rayleigh
p sig
p amp

v = GSL::Vector.linspace(0, XMAX, N)
val = GSL::Ran::rayleigh_pdf(v, sig)
val *= amp
GSL::graph(h, [v, val], "-T X -C")

h2 = GSL::Histogram.alloc(N, 0, XMAX)
for i in 0...10000 do
  x = rng.gaussian(sigma)
  y = rng.gaussian(sigma)
  r = Math.sqrt(x*x + y*y)
  h2.increment(r)
end
sig2, amp2 = h2.fit_rayleigh
p sig2
p amp2
val2 = GSL::Ran::rayleigh_pdf(v, sig2)
val2 *= amp2
GSL::graph(h2, [v, val2], "-T X -C")


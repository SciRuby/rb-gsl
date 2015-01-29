#!/usr/bin/env ruby
require("gsl")

N = 5000
BINS = 50
rng = GSL::Rng.alloc(1)

h = GSL::Histogram3d.alloc(BINS, [-5, 5], BINS, [-4, 8], BINS, [-8, 3])
sig1 = 1
sig2 = 2
sig3 = 1.3

for i in 0...N do
  r1 = rng.gaussian(sig1) + 1
  r2 = rng.gaussian(sig2) + 1.5
  r3 = rng.gaussian(sig3) - 2
  h.increment(r1, r2, r3)
end

hxy = h.xyproject
h1 = hxy.xproject
h2 = hxy.yproject
hxz = h.xzproject
h3 = hxz.yproject

x = GSL::Vector.linspace(-7, 7, 100)
a = h1.fit_gaussian   # a[0]: sigma, a[1]: mean, a[2]: height
y1 = a[2]*GSL::Ran::gaussian_pdf(x-a[1], a[0])
a = h2.fit_gaussian   # a[0]: sigma, a[1]: mean, a[2]: height
y2 = a[2]*GSL::Ran::gaussian_pdf(x-a[1], a[0])
a = h3.fit_gaussian   # a[0]: sigma, a[1]: mean, a[2]: height
y3 = a[2]*GSL::Ran::gaussian_pdf(x-a[1], a[0])

GSL::graph(h1, h2, h3, [x, y1], [x, y2], [x, y3])

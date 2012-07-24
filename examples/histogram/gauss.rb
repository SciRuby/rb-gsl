#!/usr/bin/env ruby
require("gsl")

N = 10000
MAX = 8
rng = GSL::Rng.alloc(2)

data = GSL::Ran.gaussian(rng, 1.5, N) + 2
h = GSL::Histogram.alloc(100, [-MAX, MAX])
h.increment(data)

sigma, mean, height, = h.fit_gaussian

x = GSL::Vector.linspace(-MAX, MAX, 100)
y = height*GSL::Ran::gaussian_pdf(x-mean, sigma)
GSL::graph(h, [x, y], "-T X -C -g 3")

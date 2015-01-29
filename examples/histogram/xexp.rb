#!/usr/bin/env ruby
require("gsl")

N = 10000
rng = GSL::Rng.alloc
data = rng.gamma(2, 1.5, N)
h = GSL::Histogram.alloc(100, [0, 15])
h.fill(data)
p h.bin

#result = h.fit_exponential
result = h.fit("xexp")
b = result[0]
amp = result[1]
p amp
p 1.0/b
x = GSL::Vector.linspace(0, 15, 100)
y = amp*x*GSL::Sf::exp(-x*b)
GSL::graph(h, [x, y], "-C -g 3")



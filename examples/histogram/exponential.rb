#!/usr/bin/env ruby
require("gsl")

N = 10000
rng = GSL::Rng.alloc
data = rng.exponential(2, N)
h = GSL::Histogram.alloc(100, [0, 15])
h.fill(data)

#result = h.fit_exponential
result = h.fit("exponential")
a = result[0]
b = result[1]

x = GSL::Vector.linspace(0, 15, 100)
y = a*GSL::Sf::exp(x*b)
GSL::graph(h, [x, y], "-C -g 3")



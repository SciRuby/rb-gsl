#!/usr/bin/env ruby
require("gsl")
GSL::Rng.env_setup()

N = 10000
rng = GSL::Rng.alloc

h = GSL::Histogram.alloc(20, [0, 20])

=begin
for i in 0...N do
  r = rng.poisson(5)
#  r = GSL::Ran::poisson(rng, 5)
  h.increment(r)
end
=end

v = rng.poisson(5, N)
h.fill(v)

h.normalize!
x = GSL::Vector.linspace(0, 20, 100)
y = GSL::Ran::poisson_pdf(x, 5)
GSL::graph(h, [x, y], "-C -g 3 -L 'Poisson distribution, mu = 5'")




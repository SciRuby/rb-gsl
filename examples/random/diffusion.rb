#!/usr/bin/env ruby
# 1-dimensional random work:
# This demonstrates M = 1000 trials of N = 6, 50 and 100 left-or-right steps.
# The distribution of the end points of the trials
# will be Gaussian of standard deviation sqrt(N).

require("gsl")

M = 1000
GSL::Rng.env_setup()
T = GSL::Rng::DEFAULT
seed = 2
rng = GSL::Rng.alloc(T, seed)

h = Array.new(3)
h[0] = GSL::Histogram.alloc(61, -30, 30)
h[1] = GSL::Histogram.alloc(61, -30, 30)
h[2] = GSL::Histogram.alloc(61, -30, 30)

i = 0
for n in [6, 50, 100] do
  M.times do
    s = 0
    n.times do
      ds = rng.get%2 == 0 ? 1 : -1
      s += ds
    end
    h[i].increment(s)
  end
  i += 1
end

#GSL::graph(h[0].shift(250), h[1].shift(100), h[2])
GSL::graph(h[0] + 250, h[1] + 100, h[2])

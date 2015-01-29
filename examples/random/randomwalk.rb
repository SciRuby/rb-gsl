#!/usr/bin/env ruby
# 1-dimensional random work:
# This demonstrates 1000 trials of 50 left-or-right steps.
# The distribution of the end points of the trials
# will be Gaussian of standard deviation sqrt(50).

require("gsl")

N = 50
M = 1000
GSL::Rng.env_setup()
T = GSL::Rng::DEFAULT
seed = 2
rng = GSL::Rng.alloc(T, seed)

sigma = Math::sqrt(N).to_i

h = GSL::Histogram.alloc(8*sigma+1, [-4*sigma-0.5, 4*sigma+0.5])

M.times do
  s = 0
  N.times do
    ds = rng.get%2 == 0 ? 1 : -1
    s += ds
  end
  h.increment(s)
end

x = GSL::Vector.linspace(-40, 40, 80)
y = GSL::Ran::gaussian_pdf(x, sigma)*M*2
# Factor 2 is not important, but necessary
# because only the even ranges are filled:
#   a + b = N     a: positive steps, b: negative steps
#   a - b = s     s: the end point after the N steps
# Since N = 1000 and a, b, s are integer, s must be even.
GSL::graph(h, [x, y], "-C -x #{-4*sigma} #{4*sigma}")


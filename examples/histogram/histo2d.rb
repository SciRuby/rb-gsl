#!/usr/bin/env ruby
require("gsl")
N = 10000
BINS = 100

rng = GSL::Rng.alloc("mt19937", 123456)

h2 = GSL::Histogram2d.alloc(BINS, [-8, 8], BINS, [-8, 8])

sig1 = 0.8
sig2 = 2.0

for i in 0...N do
  r1 = rng.gaussian(sig1) + 2.5
  r2 = rng.gaussian(sig2) - 1
  h2.increment(r1, r2)
end

hx = h2.xproject
hy = h2.yproject
printf("%f %f %f %f\n", h2.xmean, h2.ymean, hx.mean, hy.mean)
printf("%f %f %f %f\n", h2.xsigma, h2.ysigma, hx.sigma, hy.sigma)

x = GSL::Vector.linspace(-8, 8, 100)
result = hx.fit_gaussian
y1 = result[2]*GSL::Ran::gaussian_pdf(x-result[1], result[0])
result = hy.fit_gaussian
y2 = result[2]*GSL::Ran::gaussian_pdf(x-result[1], result[0])
GSL::graph(hx, hy, [x, y1], [x, y2], "-T X -C -g 3")



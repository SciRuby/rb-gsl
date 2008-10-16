#!/usr/bin/env ruby
require("gsl")

N = 100000
rng = GSL::Rng.alloc

v = rng.gaussian(2, N)

p v.mean
p v.variance
p v.sd

p v.variance_m(0)
p v.sd_m(0)

p v.variance_with_fixed_mean(0)
p v.sd_with_fixed_mean(0)

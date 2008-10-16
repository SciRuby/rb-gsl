#!/usr/bin/env ruby
require("gsl")

N = 10000
r = GSL::Rng.alloc
na = r.gaussian(1.0, N).to_na    # Generate N random numbers
p na.class
p na.rank
p na.size
p na.min
p na.max
h = GSL::Histogram.alloc(50, [-4, 4])
h.fill(na)
h.graph("-T X -C -g 3")

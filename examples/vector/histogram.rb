#!/usr/bin/env ruby
require("gsl")

N = 10000
r = GSL::Rng.alloc
v = r.gaussian(1.0, N)    # Generate N random numbers
h = v.histogram(50, [-4, 4])
h.graph("-T X -C -g 3")



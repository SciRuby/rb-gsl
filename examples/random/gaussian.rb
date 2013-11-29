#!/usr/bin/env ruby
require("gsl")

N = 10000
r = GSL::Rng.alloc
v = r.gaussian(1, N)
h = v.histogram(100, -4, 4)
h.graph


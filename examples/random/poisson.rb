#!/usr/bin/env ruby
require("gsl")

N = 10000
r = GSL::Rng.alloc
v = r.poisson(5, N)
h = v.histogram(20, 0, 20)
h.graph


#!/usr/bin/env ruby
require("gsl")

r = GSL::Rng.alloc
m = GSL::Matrix.rand(3, 3, r)
p m
m = GSL::Matrix.randn(3, 3, r)
p m

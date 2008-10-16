#!/usr/bin/env ruby
require("gsl")
puts("v = Vector.logspace2(1, 10000, 10)")
v = GSL::Vector.logspace2(1, 10000, 10)
p v
v.graph("-C -Y v -l y -S 4 -L 'Vector.logspace2(1, 10000, 10)'")

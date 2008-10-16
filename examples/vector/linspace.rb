#!/usr/bin/env ruby
require("gsl")
puts("v = GSL::Vector.linspace(0, 5, 10)")
v = GSL::Vector.linspace(0, 5, 10)
p v
v.graph("-C -Y v -S 4 -L 'Vector.linspace(0, 5, 10)'")

#!/usr/bin/env ruby
require("gsl")

v = GSL::Vector.alloc([1, 2, 3, 4])
p v
v.fwrite("a.dat")

v2 = GSL::Vector.alloc([5, 6, 7, 8])
p v2
v2.fprintf("b.dat")


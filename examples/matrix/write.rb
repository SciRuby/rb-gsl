#!/usr/bin/env ruby
require("gsl")
m = GSL::Matrix.alloc([1, 2], [3, 4])
m.fwrite("a.dat")

m2 = GSL::Matrix.alloc([5, 6], [7, 8])
m2.fprintf("b.dat")

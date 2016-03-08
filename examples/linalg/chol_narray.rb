#!/usr/bin/env ruby
require("gsl")
include GSL
include Linalg

m = NArray[[4.0, 2], [2, 3]]
c = Cholesky.decomp(m)
puts "decomp ->"
p c

b = NArray[1.0, 2]
puts "solve ->"
p Cholesky.solve(c, b)    # Expected [-0.125, 0.75]

b = NArray[1.0, 2]
Cholesky.svx(c, b)
puts "svx ->"
p b

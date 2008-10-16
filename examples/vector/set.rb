#!/usr/bin/env ruby
require("gsl")

v = GSL::Vector.alloc(1, 2, 3, 4, 5)

p v
p v[3]

v.set_all(9)
p v
v.set_zero
p v

v.set_basis(3)
p v

v.set(2, 5)
p v

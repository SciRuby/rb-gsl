#!/usr/bin/env ruby
require("gsl")

a = GSL::Matrix.alloc([1, 2, 3, 4, 5, 6, 7, 8, 9], 3, 3)

a[1][2] = 99.9
p a
 
a.set_all(5)
p a

a.set_zero
p a

a.set_identity
p a

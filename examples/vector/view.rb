#!/usr/bin/env ruby
require("gsl")

v = GSL::Vector.alloc(1, 2, 3, 4, 5, 6, 7, 8, 9)
p v
vv = v.subvector(2, 3)
p vv
vv.set(1, 9)
p vv
p v

m = v.matrix_view(3, 3)
p m

v2 = v.view
p v2
p v2.class

v2 = v.view(3)
p v2
p v2.size
p v2.class

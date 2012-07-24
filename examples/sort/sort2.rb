#!/usr/bin/env ruby
require("gsl")

v = GSL::Vector.alloc(1, 2, 6, 7, 8, 9, 3, 4, 5)
p v.smallest(5)
p v.largest(3)

p v.smallest_index(3)
p v.largest_index(5)

p v.sort
p v

v.sort!
p v


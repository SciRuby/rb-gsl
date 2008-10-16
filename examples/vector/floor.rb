#!/usr/bin/env ruby
require("gsl")

v = GSL::Vector[1.1, 2.7, 3.5, 5.8, -1.2, -2.8, -3.5]
p v
puts("floor")
p v.floor
puts("ceil")
p v.ceil
puts("round")
p v.round

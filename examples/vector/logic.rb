#!/usr/bin/env ruby
require("gsl")

a = GSL::Vector[1, 0, 3, 0]
b = GSL::Vector[3, 4, 0, 0]
p a.and(b)
p a.or(b)
p a.xor(b)
#p a.not

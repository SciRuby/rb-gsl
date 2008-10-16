#!/usr/bin/env ruby
require("gsl")

v1 = GSL::Vector[1..2]
v2 = GSL::Vector[3..4]
v3 = GSL::Vector[5..7]

puts("v1 = [ #{v1.to_a.join(' ')} ]")
puts("v2 = [ #{v2.to_a.join(' ')} ]")
puts("v3 = [ #{v3.to_a.join(' ')} ]\n")

puts("Ex1: v1.connect(v2, v3)")
p v1.connect(v2, v3)

puts("Ex2: GSL::Vector.connect(v1, v2, v3)")
p GSL::Vector.connect(v1, v2, v3)


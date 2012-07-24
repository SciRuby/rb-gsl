#!/usr/bin/env ruby
require("gsl")

puts("\nGSL::Vector")
v = GSL::Vector.alloc(1..9)
p v

puts("\nGSL::Vector ---> NArray")
vna = v.to_na
p vna

puts("\nNArray ---> GSL::Vector")
v2 = GSL::Vector.to_gv(vna)
p v2

puts("\nGSL::Matrix")
m = GSL::Matrix.alloc(1..9, 3, 3)
p m

puts("\nGSL::Matrix ---> NArray")
mna = m.to_na
p mna

puts("\nNArray ---> GSL::Matrix")
m2 = GSL::Matrix.to_gm(mna)
p m2



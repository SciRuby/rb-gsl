#!/usr/bin/env ruby
require("gsl")

#z = GSL::Vector[1, 0, 0, 5, 2, 0, 9]
#z = GSL::Vector[1, 0, 0, 5, 2, 0, 9].block
z = GSL::Vector::Int[1, 0, 0, 5, 2, 0, 9]
#z = GSL::Vector::Int[1, 0, 0, 5, 2, 0, 9].block
#z = GSL::Vector[1, 0, 0, 5, 2, 0, 9].and(1)

p z

puts("where != 0")
p z.where
puts("where >= 2")
p z.where { |elm| elm >= 2 }

puts("-----")

puts("where2")
p z.where2
puts("where2,  >= 2")
p z.where2 { |elm| elm >= 2 }

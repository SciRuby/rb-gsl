#!/usr/bin/env ruby
require("gsl")

x = GSL::Vector::Int[1, 2, 3]
y = GSL::Vector::Int[1, 2, 5]
z = GSL::Vector::Int[0, 2, 9]

puts("x = #{x.to_s}")
puts("y = #{y.to_s}")
puts("z = #{z.to_s}")

puts("Test x.eq(y)")
p x.eq(y)
puts("Test x.ne(y)")
p x.ne(y)
puts("Test x.ge(y)")
p x.ge(y)
puts("Test x.lt(z)")
p x.lt(z)

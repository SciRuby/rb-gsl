#!/usr/bin/env ruby
require("gsl")

x = GSL::Vector::Int[1, 2, 3]
y = GSL::Vector::Int[1, 0, 5]
a = GSL::Vector::Int[0, 0, 0]

puts("x = #{x.to_s}")
puts("y = #{y.to_s}")
puts("a = #{a.to_s}")

puts("x.any? = #{x.any?}, x.all? = #{x.all?}, x.none? = #{x.none?}")
puts("y.any? = #{y.any?}, y.all? = #{y.all?}, y.none? = #{y.none?}")
puts("a.any? = #{a.any?}, a.all? = #{a.all?}, a.none? = #{a.none?}")

puts("x.any? { |val| val > 5 } ---> #{x.any? { |val| val > 5 }}")
puts("x.any? { |val| val > 2 } ---> #{x.any? { |val| val > 2 }}")
puts("x.all? { |val| val >= 1 } ---> #{x.all? { |val| val >= 1 }}")
puts("x.all? { |val| val >= 2 } ---> #{x.all? { |val| val >= 2 }}")
puts("x.none? { |val| val == 1 } ---> #{x.none? { |val| val == 1 }}")
puts("x.none? { |val| val == 5 } ---> #{x.none? { |val| val == 5 }}")

#!/usr/bin/env ruby
require("gsl")

a = GSL::Matrix.alloc([1, 2, 3, 4, 5, 6, 7, 8, 9], 3, 3)
b = GSL::Matrix.alloc([6, 7, 8], [2, 3, 4], [3, 4, 5])

p a * 2

p 2 * a

a *= 2
p a

a /= 2
p a

p a*b

p a.mul_elements(b)

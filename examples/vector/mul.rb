#!/usr/bin/env ruby
require("gsl")

a = GSL::Vector.alloc(1, 2, 3, 4, 5)
b = GSL::Vector.alloc(6, 7, 8, 9, 10)

p a * b
p a

a *= b
p a

a /= b
p a

p a * 2
p a

a *= 2
p a

a /= 2
p a

p 2 * a
p a

p 5 / a
p a


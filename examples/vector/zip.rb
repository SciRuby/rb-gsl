#!/usr/bin/env ruby
require("gsl")
require("gsl/gplot")

a = GSL::Vector[0..4]
p a
b = GSL::Vector[2, 3, 4]
c = GSL::Vector[5, 7, 4, 8, 9, 2]

p a.zip(b, c)

p GSL::Vector.zip(a, b, c)

aa = a.to_complex
bb = b.to_complex
cc = c.to_complex

p aa.zip(bb, cc)
p GSL::Vector::Complex.zip(aa, bb, cc)

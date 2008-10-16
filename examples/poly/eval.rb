#!/usr/bin/env ruby
require("gsl")

# f(x) = 1 + 2x + 3x^3
# Evaluate at x = 2 ---> f(x) = 17
p GSL::Poly.eval([1, 2, 3], 2)
p GSL::Poly.eval([1, 2, 3].to_gv, 2)
p GSL::Poly.eval(NArray[1.0, 2, 3], 2)

# f(1) = 6, f(2) = 17, f(3) = 34
p GSL::Poly.eval([1, 2, 3], [1, 2, 3])
p GSL::Poly.eval([1, 2, 3], [1, 2, 3].to_gv)
p GSL::Poly.eval([1, 2, 3], NArray[1.0, 2, 3])


v = GSL::Vector[1, 2, 3]
p GSL::Poly.eval(v, [1, 2, 3])
p GSL::Poly.eval(v, [1, 2, 3].to_gv)
p GSL::Poly.eval(v, NArray[1.0, 2, 3])

v = NArray[1.0, 2, 3]
p GSL::Poly.eval(v, [1, 2, 3])
p GSL::Poly.eval(v, [1, 2, 3].to_gv)
p GSL::Poly.eval(v, NArray[1.0, 2, 3])

v = GSL::Vector[1, 2, 3]
x = GSL::Matrix.alloc(1...9, 3, 3)
p GSL::Poly.eval(v, x)

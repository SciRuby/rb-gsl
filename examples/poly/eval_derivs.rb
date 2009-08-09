#!/usr/bin/env ruby
require("gsl")

poly = [1, 2, 3]
p GSL::Poly::eval_derivs(poly, 1)

poly = NArray[1.0, 2, 3]
p GSL::Poly::eval_derivs(poly, 1)

poly = GSL::Vector.alloc([1, 2, 3])
p GSL::Poly::eval_derivs(poly, 1)

p poly.eval_derivs(1, 3)


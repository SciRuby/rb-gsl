#!/usr/bin/env ruby
require("gsl")

poly = [1, 2, 3]
p GSL::Poly::eval_derivs(poly, 1)   # Returned Array

poly = NArray[1.0, 2, 3]
p GSL::Poly::eval_derivs(poly, 1)   # Returned NArray

poly = GSL::Poly.alloc([1, 2, 3])   # Returned GSL::Poly
p GSL::Poly::eval_derivs(poly, 1)

p poly.eval_derivs(1, 3)


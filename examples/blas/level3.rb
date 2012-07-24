#!/usr/bin/env ruby
require("gsl")
include GSL::Blas

a = GSL::Matrix.alloc([1, 5], [4, 6])
b = GSL::Matrix.alloc([10, 2], [3, 8])
c = GSL::Matrix.alloc(2, 2)

p GSL::Blas.dgemm(NoTrans, NoTrans, 1, a, b, 1, c)
p GSL::Blas.dgemm(a, b)
p a*b


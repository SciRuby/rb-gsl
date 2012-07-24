#!/usr/bin/env ruby
require("gsl")

include GSL::Eigen

m = GSL::Matrix::Complex.alloc(2, 2)
m.set(0, 1, [0, -1])
m.set(1, 0, [0, 1])
p m

val, vec = m.eigen_hermv
p val
p vec.unpack

m2 = GSL::Matrix.alloc([1, 0, 0], [0, 0, 1], [0, 1, 0])
val, vec = m2.eigen_symmv
p val
p vec.unpack





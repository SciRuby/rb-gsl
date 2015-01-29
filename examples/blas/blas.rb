#!/usr/bin/env ruby
require("gsl")

a = GSL::Vector[0.11, 0.12, 0.13, 0.21, 0.22, 0.23]
b = GSL::Vector[1011, 1012, 1021, 1022, 1031, 1032]

A = a.matrix_view(2, 3)
B = b.matrix_view(3, 2)

C = GSL::Blas.dgemm(GSL::Blas::NoTrans, GSL::Blas::NoTrans, 1.0, A, B, 0.0)
p C

p A*B

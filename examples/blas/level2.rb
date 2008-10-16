#!/usr/bin/env ruby
require("gsl")
include GSL::Blas

x = GSL::Vector[1, 2, 3]
y = GSL::Vector[4, 5, 6]

A = GSL::Matrix.alloc([1, 2, 3], [4, 5, 6], [7, 8, 9])
# [28, 64, 100]
p GSL::Blas.dgemv(CblasNoTrans, 2, A, x)
p dgemv(CblasNoTrans, 2, A, x)

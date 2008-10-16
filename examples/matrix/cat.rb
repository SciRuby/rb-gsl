#!/usr/bin/env ruby
require("gsl")

a = GSL::Matrix::Int[1..4, 2, 2]
b = GSL::Matrix::Int[5..10, 2, 3]
p a.horzcat(b)
p GSL::Matrix::Int.horzcat(a, b)

a = GSL::Matrix::Int[1..4, 2, 2]
b = GSL::Matrix::Int[5..10, 3, 2]
p b
p a.vertcat(b)
p GSL::Matrix::Int.vertcat(a, b)

#!/usr/bin/env ruby
require("gsl")
require("../gsl_test2.rb")
include GSL::Test
include Linalg

m = GSL::Matrix.pascal(6)

p m
c_exp = GSL::Matrix[[1, 0, 0, 0, 0, 0],
               [1, 1, 0, 0, 0, 0],
               [1, 2, 1, 0, 0, 0],
               [1, 3, 3, 1, 0, 0],
               [1, 4, 6, 4, 1, 0],
               [1, 5, 10, 10, 5, 1]]

c = m.cholesky_decomp
a = c.lower
test2(a == c_exp, "#{m.class}#cholesky_decomp") 
test2((a*a.trans) == m, "#{m.class}#cholesky_decomp")

#!/usr/bin/env ruby
require("gsl")
require("../gsl_test2.rb")
include GSL::Test
include Linalg

A = GSL::Matrix[1..4, 2, 2]
I = GSL::Matrix.identity(2)
Ainv = A.inv

u, v, s = A.SV_decomp
sm = s.to_m_diagonal
sinv = s.collect { |elm| 1.0/elm }.to_m_diagonal
a = u*sm*v.trans
ainv = v*sinv*u.trans
test2(a == A, "#{A.class}#SV_decomp")
test2(ainv == Ainv, "#{A.class}#SV_decomp")

test2((u.trans*u) == I, "#{A.class}#SV_decomp")
test2((v.trans*v) == I, "#{A.class}#SV_decomp")
test2(A*v == u*sm, "#{A.class}#SV_decomp")
test2(A.trans*u == v*sm, "#{A.class}#SV_decomp")



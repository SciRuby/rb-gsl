#!/usr/bin/env ruby
require("gsl")
require("../gsl_test.rb")
include GSL::Test

dbleps = 1e-6
n = 1
x = GSL::Vector.alloc(0.002)
y = GSL::Vector.alloc(-0.921)
expected = GSL::Vector.alloc(0.002)
Blas.dcopy(x, y)
for i in 0...n do
  GSL::Test::test_rel(y[i], expected[i], dbleps, "dcopy")
end

x = GSL::Vector::Complex.alloc([[ 0.315, -0.324]])
y = GSL::Vector::Complex.alloc([[-0.312, -0.748]])
expected = GSL::Vector::Complex.alloc([[0.315, -0.324]])
Blas.zcopy(x, y)
for i in 0...n do
  GSL::Test::test_rel(y[i].re, expected[i].re, dbleps, "zcopy real")
  GSL::Test::test_rel(y[i].im, expected[i].im, dbleps, "zcopy imag")
end

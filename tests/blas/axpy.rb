#!/usr/bin/env ruby
require("gsl")
require("../gsl_test.rb")
include GSL::Test

dbleps = 1e-6
n = 1
alpha = -0.3
x = GSL::Vector.alloc(0.029)
y = GSL::Vector.alloc(-0.992)
expected = GSL::Vector.alloc(-1.0007)
y2 = Blas.daxpy(alpha, x, y)
for i in 0...n do
  GSL::Test::test_rel(y2[i], expected[i], dbleps, "daxpy")
end

alpha = GSL::Complex.alloc(0, 1)
x = GSL::Vector::Complex.alloc([[0.776, -0.671]])
y = GSL::Vector::Complex.alloc([[0.39, 0.404]])
expected = GSL::Vector::Complex.alloc([[1.061, 1.18]])
y2 = Blas.zaxpy(alpha, x, y)
for i in 0...n do
  GSL::Test::test_rel(y2[i].re, expected[i].re, dbleps, "zaxpy real")
  GSL::Test::test_rel(y2[i].im, expected[i].im, dbleps, "zaxpy imag")
end

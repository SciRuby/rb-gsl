#!/usr/bin/env ruby
require("gsl")
require("./gsl_test.rb")
include GSL::Test
include Math

for i in 0...10
  r = (i - 5.0)*0.3
  t = 2.0*M_PI*i/5.0
  x = r*cos(t)
  y = r*sin(t)
  z = GSL::Complex.polar(r, t)
  desc = sprintf("gsl_complex_polar real part at (r=%g,t=%g)", r, t)
  GSL::Test.test_rel(z.real, x, 10*GSL::DBL_EPSILON, desc)
  desc = sprintf("gsl_complex_polar imag part at (r=%g,t=%g)", r, t)
  GSL::Test.test_rel(z.imag, y, 10*GSL::DBL_EPSILON, desc)
end

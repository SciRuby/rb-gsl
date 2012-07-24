#!/usr/bin/env ruby
require 'gsl'
include GSL
include Math

f454 = Function.alloc{ |x|
  x2 = x*x
  x3 = x2*x
  x3*log(((x2-1)*(x2-2)).abs)
}
exp_result = 5.274080611672716401E+01
exp_abserr = 1.755703848687062418E-04
pts = [0, 1, sqrt(2), 3]
result = f454.qagp(pts, 0.0, 1e-3)
p result

puts("exp_result: #{exp_result}")
puts("exp_abserr: #{exp_abserr}")

p Integration.qagp(f454, pts, [0.0, 1e-3])

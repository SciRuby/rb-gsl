#!/usr/bin/env ruby
require 'gsl'
include GSL
include Math

f459 = Function.alloc { |x|
  1.0/(5.0*x*x*x + 6.0)
}

exp_result = -8.994400695837000137E-02
exp_abserr = 1.185290176227023727E-06

result = f459.qawc([-1.0, 5.0], 0, [0.0, 1e-3])
p result
puts("exp_result: #{exp_result}")
puts("exp_abserr: #{exp_abserr}")

p Integration.qawc(f459, [-1.0, 5.0], 0)

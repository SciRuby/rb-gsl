#!/usr/bin/env ruby
require 'gsl'
include GSL
include Math

f456 = Function.alloc { |x|
  if x.zero?
    val = 0.0
  else
    val = log(x)
  end
  val
}

exp_result = -1.281368483991674190E-01
exp_abserr =  6.875028324415666248E-12

table = [10.0*PI, 1.0, GSL::Integration::SINE, 1000]
result = f456.qawo(0.0, [0.0, 1e-7], table)
p result
puts("exp_result: #{exp_result}")
puts("exp_abserr: #{exp_abserr}")

table = Integration::QAWO_Table.alloc(10.0*PI, 1.0, GSL::Integration::SINE, 1000)
p f456.qawo(0.0, [0.0, 1e-7], table)
p f456.qawo(0.0, table)

p Integration.qawo(f456, 0.0, table)
p Integration.qawo(f456, 0.0, [0.0, 1e-7], table)

#!/usr/bin/env ruby
require 'gsl'
include GSL
include Math

f458 = Function.alloc { |x|
  if x.zero?
    val = 0.0
  else
    u = log(x)
    v = 1.0 + u*u
    val = 1.0/(v*v)
  end
  val
}

exp_result = -1.892751853489401670E-01
exp_abserr =  1.129133712015747658E-08

table = [0.0, 0.0, 1, 0]
result = f458.qaws([0.0, 1.0], table, [0.0, 1e-7])
p result
puts("exp_result: #{exp_result}")
puts("exp_abserr: #{exp_abserr}")

table = Integration::QAWS_Table.alloc(0.0, 0.0, 1, 0)
result = f458.qaws([0.0, 1.0], table, [0.0, 1e-7])
p result

p Integration.qaws(f458, [0.0, 1.0], table)

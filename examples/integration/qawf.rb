#!/usr/bin/env ruby
require 'gsl'
include GSL
include Math

f457 = Function.alloc { |x|
  if x.zero?
    val = 0.0
  else
    val = 1.0/sqrt(x)
  end
  val
}

exp_result =  9.999999999279802765E-01
exp_abserr =  1.556289974669056164E-08

table = [PI/2.0, 1.0, GSL::Integration::COSINE, 1000]
result = f457.qawf(0.0, 1e-7, table)
p result
puts("exp_result: #{exp_result}")
puts("exp_abserr: #{exp_abserr}")

w = Integration::Workspace.alloc
wc = Integration::Workspace.alloc

limit = 1000
table = Integration::QAWO_Table.alloc(PI/2.0, 1.0, GSL::Integration::COSINE, 1000)
p f457.qawf(0.0, table)
p f457.qawf(0.0, 1e-7, table)
p f457.qawf(0.0, 1e-7, limit, table)
p f457.qawf(0.0, limit, table)
p f457.qawf(0.0, 1e-7, limit, w, wc, table)
p f457.qawf(0.0, w, wc, table)
p f457.qawf(0.0, limit, w, wc, table)
#p f457.qawf(0.0, limit, w, table)

p Integration.qawf(f457, 0.0, table)
p Integration.qawf(f457, 0.0, 1e-7, table)
p Integration.qawf(f457, 0.0, 1e-7, limit, table)
p Integration.qawf(f457, 0.0, limit, w, wc, table)

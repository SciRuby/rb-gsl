#!/usr/bin/env ruby
require 'gsl'
include GSL
include Math

f = Function.alloc{ |x|
  exp(-x*x)
}

exact = sqrt(PI)

result, = f.qagi
puts("QAGI")
puts("exp(-x*x), x = -infty --- +infty")
printf("exact  = %.18f\n", exact)
printf("result = %.18f\n\n", result)

p Integration.qagi(f)

w = Integration::Workspace.alloc(1000)
xmin = 0.0

puts("QAGIU")
result, = f.integration_qagiu(xmin, [0, 1e-6], w)
puts("exp(-x*x), x = 0 --- +infty")
printf("exact  = %.18f\n", exact/2)
printf("result = %.18f\n", result)
p w.to_a

p Integration.qagiu(f, xmin)

puts("QAGIL")
result, = f.integration_qagil(0.0, 0.0, 1e-7, 1000, w)
puts("exp(-x*x), x = -infty --- 0")
printf("exact  = %.18f\n", exact/2)
printf("result = %.18f\n\n", result)

f455 = Function.alloc { |x|
  log(x)/(1.0 + 100.0*x*x)
}

exp_result = -3.616892186127022568E-01
exp_abserr = 3.016716913328831851E-06

result = f455.qagiu(0.0, [0.0, 1e-3])
p result
puts("exp_result: #{exp_result}")
puts("exp_abserr: #{exp_abserr}")


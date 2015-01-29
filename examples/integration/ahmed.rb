#!/usr/bin/env ruby
# Ahmed's Integral
# See e.g. http://mathworld.wolfram.com/AhmedsIntegral.html
#
require("gsl")
include GSL
include Math

f = Function.alloc { |x|
  sqrtx22 = sqrt(x*x + 2)
  atan(sqrtx22)/(sqrtx22*(x*x + 1))
}

val = f.qng(0, 1)[0]

puts("Ahmed's integral")
puts("Expect: 5pi^2/96 = #{5.0*M_PI*M_PI/96}")
puts("QNG result:        #{val}")




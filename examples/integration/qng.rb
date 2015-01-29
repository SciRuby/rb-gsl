#!/usr/bin/env ruby
require("gsl")
include GSL
include Math

f = Function.alloc { |x| Math::sin(x)/x }

p f.qng(0, 2.0*Math::PI)

# Singular at x = 0
f2 = Function.alloc { |x| exp(-x)/sqrt(x) }

# This will fail...
#p f2.qng(0, 1)

p Integration.qng(f, [0, 2*Math::PI])


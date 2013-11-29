#!/usr/bin/env ruby
require("gsl")
include GSL

# Singular at x = 0
f = Function.alloc { |x| 1.0/Math::sqrt(x) }

puts("QAGS")
p f.qags(0, 1)
p f.qags([0, 1])

p Integration.qags(f, 0, 1)
p Integration.qags(f, [0, 1])


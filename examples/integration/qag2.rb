#!/usr/bin/env ruby
require("gsl")

# Singular at x = 0

f = GSL::Function.alloc { |x| 1.0/Math::sqrt(x) }

puts("QAG")
p f.qag(0, 1)
p f.qag([0, 1])

p GSL::Integration.qag(f, 0, 1)
p GSL::Integration.qag(f, [0, 1])


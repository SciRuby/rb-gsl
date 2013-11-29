#!/usr/bin/env ruby
require 'gsl'

a = GSL::Complex.rect(1, 2)

p 2 + a
p 2 - a
p 2*a

p a - 3
p a * 3

p 2.0/a  # 0.4 - 0.8i


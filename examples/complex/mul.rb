#!/usr/bin/env ruby
require 'gsl'

a = GSL::Complex.rect(1, 2)
b = GSL::Complex[3, 4]

p a*b  # -5 + 10i
p a.mul(b)
p a

p a.div(b)
p a

a *= b
p a

a /= b
p a

p a.mul(2)
p a.mul_real(2)
p a.div_real(2)

p a.conjugate
p a.inverse
p 1.0/a
p a.negative
p -a

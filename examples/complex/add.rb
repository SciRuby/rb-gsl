#!/usr/bin/env ruby
require 'gsl'

a = GSL::Complex.rect(1, 2)
b = GSL::Complex[3, 4]

p a + b
p a

a += b
p a

p a - b
p a

a -= b
p a

p a + b + a

p a.add(b)
p a

p a.sub(b)
p a

p a + 5
p a.add_real(5)
p a

p a.sub_real(5)
p a

p a.add_real(6).sub_real(6)



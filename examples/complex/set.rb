#!/usr/bin/env ruby
require 'gsl'

a = GSL::Complex.rect(1, 2)
p a

a.set(5)
p a

a.set(3, 6)
p a

p a.re
p a.im
p a.real
p a.imag
p a.REAL
p a.IMAG

a.re = 7
p a

a.im = 1
p a

a.real = 4
p a

a.imag = 9
p a

a.set_real(2)
p a

a.set_imag(3)
p a


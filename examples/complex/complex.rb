#!/usr/bin/env ruby
require 'gsl'

a = GSL::Complex.rect(1, 2)
b = GSL::Complex[3, 4]
c = GSL::Complex.alloc(5, 6)
d = GSL::Complex.alloc([7, 8])
p a
p b
p c
p d

e = GSL::Complex.polar(1, Math::PI/6)
p e

p a.abs
p a.abs2
p a.logabs
p Math::log(a.abs)

p e.abs
p e.arg
p Math::PI/6
exit


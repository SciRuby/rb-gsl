#!/usr/bin/env ruby
require("gsl")

p Jac::jacobi_P1(0.2, 0.2, 0.5)

x = GSL::Vector.linspace(-1, 1, 6)
p Jac::jacobi(x, 1, 0.2, 0.5)

p Jac::jacobi_P1(x, 0.2, 0.5)

p Jac::jacobi_zeros(6, 0.2, 0.5)

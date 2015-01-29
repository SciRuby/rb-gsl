#!/usr/bin/env ruby
# Find synchrotron spectrum peak
require("gsl")

# Create Function object from the module function
F = GSL::Function.alloc { |x| GSL::Sf::synchrotron_1(x) }
# Derivative of the function
DF = GSL::Function.alloc { |x|
  result, abserr, status = F.deriv_central(x, 1e-6)
  result
}
# Find zero-point of the derivative in interval (0.01, 5)
peak, = DF.fsolve(0.01, 5)
printf("A peak is found at %3.2f.\n", peak)

x = GSL::Vector.linspace(0, 5, 100)
s = GSL::Sf::synchrotron_1(x)
s.graph(x, "-C -g 3 -X x -L 'Sf::synchrotron_1(x)'")

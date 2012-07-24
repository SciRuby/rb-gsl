#!/usr/bin/env ruby
require("gsl")

f = GSL::Function.alloc { |x| Math::exp(-x) }

puts("f(x) = exp(-x)");
puts("Derivative: -exp(-x)")

h = 1e-8
x = GSL::Vector.linspace(0, 5, 50)
y = f.eval(x)
dy, derr = f.deriv_central(x, h)
GSL::graph(x, y, dy, "-T X -C -g 3 -X x -L 'f(x) = exp(-x), and its derivative'")

#!/usr/bin/env ruby
require("gsl")
n = 50
x = GSL::Vector.linspace(0, 3, n)
y1 = GSL::Ran::exponential_pdf(x, 1)
y2 = GSL::Ran::exponential_pdf(x, 2)
GSL::Vector.graph(x, y1, y2, "-T X -C -g 3 -X x -Y 'p(x)' -y 0 1 -L 'Exponential distribution' --toggle-rotate-y-label")

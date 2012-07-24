#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector.linspace(0, 5, 100)
y1 = GSL::Ran::flat_pdf(x, 0.5, 2.5)
y2 = GSL::Ran::flat_pdf(x, 1.2, 4.8)
GSL::Vector.graph(x, y1, y2, "-T X -C -g 3  -X x -Y 'p(x)' -y 0 1 -L 'Flat distribution, a = 0.5, b = 2.5' --toggle-rotate-y-label")


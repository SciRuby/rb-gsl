#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector.linspace(0, 5, 100)
y1 = GSL::Ran::gamma_pdf(x, 1, 1)
y2 = GSL::Ran::gamma_pdf(x, 2, 1)
y3 = GSL::Ran::gamma_pdf(x, 3, 1)
GSL::Vector.graph(x, y1, y2, y3, "-T X -C -g 3 -X x -Y 'p(x)' -L 'Gamma distribution, a = 1, 2, 3' --toggle-rotate-y-label")


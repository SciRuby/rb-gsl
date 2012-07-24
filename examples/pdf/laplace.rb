#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector.linspace(-5, 5, 100)
y1 = GSL::Ran::laplace_pdf(x, 1)
y2 = GSL::Ran::laplace_pdf(x, 2)
GSL::Vector.graph(x, y1, y2, "-T X -g 3 -C -X x -Y 'p(x)' -x -5 5 -L 'Laplace distribution, a = 1, 2' --toggle-rotate-y-label")


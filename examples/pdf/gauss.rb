#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector.linspace(-5, 5, 100)
y1 = GSL::Ran::gaussian_pdf(x, 1)
y2 = GSL::Ran::gaussian_pdf(x, 2)
GSL::Vector.graph(x, y1, y2, "-T X -C -g 3 -X x -Y 'p(x)' -x -5 5 -y 0 0.5 -L 'Gaussian distribution, sigma = 1, 2' --toggle-rotate-y-label")

#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector.linspace(0, 5, 100)
y1 = GSL::Ran::gaussian_tail_pdf(x, 1.5, 1)
y1.graph(x, "-C -X x -Y 'p(x)' -x 0 5 -L 'Gaussian-tail distribution, s = 1, a = 1.5' --toggle-rotate-y-label")

#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector.linspace(0, 5, 100)
y1 = GSL::Ran::rayleigh_tail_pdf(x, 1, 1)
y2 = GSL::Ran::rayleigh_tail_pdf(x, 0.5, 2)
GSL::Vector.graph(x, y1, y2, "-T X -g 3 -C -X x -Y 'p(x)' -x 0 5 -y 0 1.2 -L 'Rayleigh-tail distribution, a = 1, s = 1' --toggle-rotate-y-label")

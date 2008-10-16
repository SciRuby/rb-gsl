#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector.linspace(0, 3, 100)
y1 = GSL::Ran::lognormal_pdf(x, 0, 1)
y2 = GSL::Ran::lognormal_pdf(x, 1, 1)
GSL::graph(x, y1, y2, "-T X -g 3 -C -X x -Y 'p(x)' -L 'Lognormal distribution, zeta=0,1, s=1' --toggle-rotate-y-label")

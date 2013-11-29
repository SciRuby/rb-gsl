#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector.linspace(0, 1, 100)
y1 = GSL::Ran::beta_pdf(x, 2, 2)
y2 = GSL::Ran::beta_pdf(x, 4, 1)
y3 = GSL::Ran::beta_pdf(x, 1, 4)
GSL::graph(x, y1, y2, y3, "-T X -g 3 -C -X x -Y 'p(x)' -L 'Beta distribution' --toggle-rotate-y-label")

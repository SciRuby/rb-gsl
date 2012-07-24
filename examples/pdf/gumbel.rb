#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector.linspace(-2, 2, 100)
y1 = GSL::Ran::gumbel1_pdf(x, 1, 1)
y2 = GSL::Ran::gumbel2_pdf(x, 1, 1)
GSL::graph(x, y1, y2, "-T X -g 3 -C -X x -Y 'p(x)' -L 'Type 1, 2 Gumbel distribution' --toggle-rotate-y-label")

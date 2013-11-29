#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector.linspace(-4, 4, 100)
y1 = GSL::Ran::tdist_pdf(x, 1)
y2 = GSL::Ran::tdist_pdf(x, 5)
GSL::graph(x, y1, y2, "-T X -C -g 3 -X x -Y 'p(x)' -L 'Student t distribution, n=1' --toggle-rotate-y-label")

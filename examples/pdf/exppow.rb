#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector.linspace(-5, 5, 100)
y1 = GSL::Ran::exppow_pdf(x, 1, 2.5)
y2 = GSL::Ran::exppow_pdf(x, 1, 0.5)
GSL::Vector.graph(x, y1, y2, "-T X -g 3 -C -X x -Y 'p(x)' -x -5 5 -L 'Exppow distribution, a = 1, b = 2.5, 0' --toggle-rotate-y-label")

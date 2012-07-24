#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector.linspace(-5, 5, 100)
y1 = GSL::Ran::logistic_pdf(x, 1)
y2 = GSL::Ran::logistic_pdf(x, 2)
GSL::graph(x, y1, y2, "-T X -g 3 -C -X x -Y 'p(x)' -x -5 5 -y 0 0.3 -L 'Logistic distribution' --toggle-rotate-y-label")

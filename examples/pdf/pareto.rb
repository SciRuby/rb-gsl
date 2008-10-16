#!/usr/bin/env ruby
require("gsl")
n = 100
x = GSL::Vector.linspace(0, 5, n)
y1 = GSL::Ran::pareto_pdf(x, 1, 1)
y2 = GSL::Ran::pareto_pdf(x, 3, 2)
GSL::graph(x, y1, y2, "-T X -g 3 -C -x 0 5 -y 0 2 -X x -Y 'p(x)' -L 'Pareto distribution' --toggle-rotate-y-label")

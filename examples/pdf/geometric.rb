#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector[0..5]
y = GSL::Ran::geometric_pdf(x, 0.5)
y.graph_step(x, "-C -X x -Y 'p(x)' -y 0 0.7 -L 'Geometric distribution, p = 0.5' --toggle-rotate-y-label")

#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector[0..10]
y = GSL::Ran::logarithmic_pdf(x, 0.7)
y.graph_step(x, "-C -X x -Y 'p(x)' -y 0 0.7 -L 'Logarithmic distribution, p = 0.7' --toggle-rotate-y-label")

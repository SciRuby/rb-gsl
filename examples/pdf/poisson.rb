#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector[0..10]
y = GSL::Ran::poisson_pdf(x, 2.5)
y.graph_step(x, "-C -X x -Y 'p(x)' -L 'Poisson distribution, mu = 2.5' --toggle-rotate-y-label")

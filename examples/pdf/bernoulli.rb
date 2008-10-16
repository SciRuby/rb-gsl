#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector[0..1]
y = GSL::Ran::bernoulli_pdf(x, 0.7)
y.graph_step(x, "-C -X x -Y 'p(x)' -y 0 1 -L 'Bernoulli distribution, p = 0.7' --toggle-rotate-y-label")

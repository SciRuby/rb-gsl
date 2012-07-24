#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector.linspace(-5, 10, 100)
y1 = GSL::Ran::landau_pdf(x)
y1.graph(x, "-C -X x -Y 'p(x)' -L 'Landau distribution' --toggle-rotate-y-label")

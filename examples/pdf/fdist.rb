#!/usr/bin/env ruby
require("gsl")
n = 100
x = GSL::Vector.linspace(0.01, 6, n)
y1 = GSL::Ran::fdist_pdf(x, 10, 20)
y2 = GSL:Ran::fdist_pdf(x, 3, 4)
GSL::graph(x, y1, y2, "-T X -g 3 -C -y 0 1 -X x -Y 'p(x)' -L 'F distribution' --toggle-rotate-y-label")

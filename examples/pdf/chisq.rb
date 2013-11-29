#!/usr/bin/env ruby
require("gsl")
n = 50
x = GSL::Vector.linspace(0.01, 3, n)
y1 = GSL::Ran::chisq_pdf(x, 1)
y2 = GSL::Ran::chisq_pdf(x, 2)
y3 = GSL::Ran::chisq_pdf(x, 3)
GSL::graph(x, y1, y2, y3, "-T X -g 3 -C -y 0 1 -X x -Y 'p(x)' -L 'Chi^2 distribution, dof=1,2,3' --toggle-rotate-y-label")

#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector.linspace(0, 2, 100)
y1 = GSL::Ran::weibull_pdf(x, 1, 1)
y2 = GSL::Ran::weibull_pdf(x, 1, 2)
y3 = GSL::Ran::weibull_pdf(x, 2, 3)
GSL::graph(x, y1, y2, y3,  "-T X -g 3 -C -x 0 2 -y 0 1.5 -X x -Y 'p(x)' -L 'Weibull distribution' --toggle-rotate-y-label")


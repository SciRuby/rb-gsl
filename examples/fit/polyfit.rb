#!/usr/bin/env ruby
require("gsl")

x = GSL::Vector[0.0, 0.3, 0.8, 1.1, 1.6, 2.3]
y = GSL::Vector[0.5, 0.82, 1.14, 1.25, 1.35, 1.40]
coef, cov, chisq, status = GSL::MultiFit.polyfit(x, y, 2)
p coef
x2 = GSL::Vector.linspace(0, 2.5, 20)
GSL::graph([x, y], [x2, coef.eval(x2)], "-C -g 3 -S 4 -X X -Y Y")

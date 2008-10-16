#!/usr/bin/env ruby
require("gsl")
include GSL

x = Vector[1, 2, 3, 4, 5]
y = Vector[5.5, 43.1, 128, 290.7, 498.4]
coef, cov, chisq, status = Poly.fit(x, y, 3)
p coef
x2 = Vector.linspace(1, 5, 20)
graph([x, y], [x2, coef.eval(x2)], "-C -g 3 -S 4 -X X -Y Y")

#

x = Vector[0, 0.3, 0.8, 1.1, 1.6, 2.3]
y = Vector[0.5, 0.82, 1.14, 1.25, 1.35, 1.40]
coef, cov, chisq, status = MultiFit.polyfit(x, y, 2)
p coef
x2 = Vector.linspace(0, 2.5, 20)
graph([x, y], [x2, coef.eval(x2)], "-C -g 3 -S 4 -X X -Y Y")



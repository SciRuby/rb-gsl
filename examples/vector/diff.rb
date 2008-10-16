#!/usr/bin/env ruby
require("gsl")

x = GSL::Vector.linspace(1, 5, 5)
y = GSL::pow_5(x)
y1 = y.diff
y2 = y1.diff
y3 = y2.diff
y4 = y3.diff
p y.diff(2) == y2
p y.diff(3) == y3
p y.diff(4) == y4

GSL::graph(x, y, y1, y2, y3, "-C -g 3 -l x -l y -x 1 5")

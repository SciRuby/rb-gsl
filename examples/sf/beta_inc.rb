#!/usr/bin/env ruby
require("gsl")
n = 50
x = GSL::Vector.linspace(0.01, 1, n)
y1 = GSL::Sf::beta_inc(0.5, 5.0, x)
y2 = GSL::Sf::beta_inc(1.0, 3.0, x)
y3 = GSL::Sf::beta_inc(8.0, 10.0, x)
y4 = GSL::Sf::beta_inc(5.0, 0.5, x)
GSL::Vector.graph([x, y1], [x, y2], [x, y3], [x, y4], "-T X -C -g 3 -y 0 1.1 -X x -L 'Incomplete beta functions'")

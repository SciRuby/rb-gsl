#!/usr/bin/env ruby
require("gsl")
n = 100
x = GSL::Vector.linspace(0, 5, n)
y1 = GSL::Sf::gegenpoly_1(0, x)
y2 = GSL::Sf::gegenpoly_2(0, x)
y3 = GSL::Sf::gegenpoly_3(0, x)
GSL::Vector.graph(x, y1, y2, y3, "-T X -C -g 3 -X x -L 'Gegenbauer functions'")

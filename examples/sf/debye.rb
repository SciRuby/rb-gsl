#!/usr/bin/env ruby
require("gsl")
n = 50
x = GSL::Vector.linspace(0, 20, n)
y1 = GSL::Sf::debye_1(x)
y2 = GSL::Sf::debye_2(x)
y3 = GSL::Sf::debye_3(x)
y4 = GSL::Sf::debye_4(x)
GSL::Vector.graph(x, y1, y2, y3, y4, "-T X -C -g 3 -X x -L 'Debye functions'")

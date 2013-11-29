#!/usr/bin/env ruby
require("gsl")
n = 50
x = GSL::Vector.linspace(0, 2, n)

y2 = GSL::Sf::transport_2(x)
y3 = GSL::Sf::transport_3(x)
y4 = GSL::Sf::transport_4(x)
y5 = GSL::Sf::transport_5(x)
GSL::Vector.graph(x, y2, y3, y4, y5, "-T X -C -g 3 -X x -Y 'J(x)' --toggle-rotate-y-label -L 'Transport functions'")

#!/usr/bin/env ruby
require("gsl")
n = 50
x = GSL::Vector.linspace(-1, 1, n)
y1 = GSL::Sf::legendre_P1(x)
y2 = GSL::Sf::legendre_P2(x)
y3 = GSL::Sf::legendre_P3(x)
y4 = GSL::Sf::legendre_Pl(4, x)
y5 = GSL::Sf::legendre_Pl(5, x)
GSL::Vector.graph(x, y1, y2, y3, y4, y5, "-T X -C -g 3 -X x -Y 'P(x)' --toggle-rotate-y-label -L 'Legendre P functions'")

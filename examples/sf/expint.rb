#!/usr/bin/env ruby
require("gsl")
n = 100
x = GSL::Vector.linspace(0.01, 2, n)
y1 = GSL::Sf::expint_E1(x)
y2 = GSL::Sf::expint_E2(x)
yi = GSL::Sf::expint_Ei(x)
GSL::Vector.graph(x, y1, y2, yi, "-T X -C -g 3 -X x -Y 'E(x)' --toggle-rotate-y-label -L 'Exponential integral E1, E2, Ei'")

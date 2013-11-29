#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector.linspace(0.01, 12, 100)
s = GSL::Sf::Shi(x)   # s and c are GSL::Vector
c = GSL::Sf::Chi(x)
GSL::Vector.graph(x, s, c, "-T X -C -g 3 -y -1.5 2 -L 'Hyperbolic integrals Shi(x), Chi(x)'")

#!/usr/bin/env ruby
require("gsl")
n = 50
x = GSL::Vector.linspace(-0.35, 10, n)
GSL::Sf::lambert_W0(x).graph(x, "-T X -C -g 3 -X x -L 'Lambert W function'")

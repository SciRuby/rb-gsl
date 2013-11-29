#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector.linspace(0.01, 14, 50)
y = GSL::Sf::lngamma(x)
GSL::graph(x, y, "-T X -C -g 3 -l y -L 'lngamma(x)'")

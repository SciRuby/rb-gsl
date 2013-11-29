#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector.linspace(0, 20, 100)
y = GSL::Sf::dawson(x)
y.graph(x, "-T X -C -g 3 -L 'Dawson integral'")

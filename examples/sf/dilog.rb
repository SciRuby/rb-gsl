#!/usr/bin/env ruby
require("gsl")
n = 100
x = GSL::Vector.linspace(0, 20, n)
y = GSL::Sf::dilog(x)
y.graph(x, "-T X -C -g 3 -X x -L 'Dilogarithm'")

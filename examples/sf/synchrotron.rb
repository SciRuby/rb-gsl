#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector.linspace(0, 5, 100)
s = GSL::Sf::synchrotron_1(x)
s.graph(x, "-C -g 3 -X x -L 'Sf::synchrotron_1(x)'")

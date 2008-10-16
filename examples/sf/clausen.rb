#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector.linspace(0, 20, 100)
a = GSL::Sf::clausen(x)
GSL::graph(x, a, "-C -g 3 -L 'Sf::clausen(x)'")


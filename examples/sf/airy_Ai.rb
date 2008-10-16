#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector.linspace(-10, 10, 100)
a = GSL::Sf::airy_Ai(x)
a.graph(x, "-C -g 3 -L 'GSL::Sf::airy_Ai(x)'")




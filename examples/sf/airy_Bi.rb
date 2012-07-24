#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector.linspace(-1, 4, 50)
a = GSL::Sf::airy_Bi(x)
a.graph(x, "-C -g 3 -L 'Sf::airy_Bi(x)'")




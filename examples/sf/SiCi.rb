#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector.linspace(0.1, 12, 100)
s = GSL::Sf::Si(x) - Math::PI/2.0
c = GSL::Sf::Ci(x)
GSL::Vector.graph(x, s, c, "-T X -C -g 3 -y -1.5 1 -L 'Red: Sine integral, Green: Cosine integral'")

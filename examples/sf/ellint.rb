#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector.linspace(0, 0.9, 100)
k = GSL::Sf::ellint_Kcomp(x)
e = GSL::Sf::ellint_Ecomp(x)
GSL::Vector.graph(x, k, e, "-T X -g 3 -C -L 'Red: Kcomp(x), Green: Ecomp(x)'")

#!/usr/bin/env ruby
require("gsl")
s = GSL::Vector.linspace(1.05, 13, 100)
y = GSL::Sf::zetam1(s)
y.graph(s, "-C -l y -L 'Zeta function z(s)-1' -X s -Y 'z(s)-1' --toggle-rotate-y-label")

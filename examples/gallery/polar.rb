#!/usr/bin/env ruby
require("gsl")

phi = GSL::Vector.linspace(0, 2*Math::PI, 200)
r1 = GSL::Sf::sin(4*phi)
x1 = r1*GSL::Sf::cos(phi)
y1 = r1*GSL::Sf::sin(phi)
r2 = GSL::Sf::cos(4*phi)
x2 = r2*GSL::Sf::cos(phi)
y2 = r2*GSL::Sf::sin(phi)
GSL::Vector.graph([x1, y1], [x2, y2])

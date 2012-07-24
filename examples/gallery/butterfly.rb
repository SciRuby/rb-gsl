#!/usr/bin/env ruby
require("gsl")
phi = GSL::Vector.linspace(0, 12*Math::PI, 800)
r = GSL::Sf::exp(GSL::Sf::cos(phi)) - 2*GSL::Sf::cos(4*phi) + GSL::pow_5(GSL::Sf::sin(phi/12))
x = r*GSL::Sf::cos(phi)
y = r*GSL::Sf::sin(phi)
GSL::graph(x, y)

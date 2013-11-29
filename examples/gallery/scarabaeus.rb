#!/usr/bin/env ruby
require("gsl")

t = GSL::Vector.linspace(0, 2*Math::PI, 100)

a = 2
b = 3

cost = GSL::Sf::cos(t)
r = b*GSL::Sf::cos(2*t) - a*cost

x = r*cost
y = r*GSL::Sf::sin(t)
GSL::graph(x, y, "-T X -C")
